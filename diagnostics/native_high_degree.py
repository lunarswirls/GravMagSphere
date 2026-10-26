"""validate the native spectral path without shtools and stop at measured resource limits"""

import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
from time import perf_counter

if sys.argv[1:] in (['--help'],['-h']):
    print(f'usage: python {Path(__file__).name} [degree]\n'
          'validate one degree or run the native high-degree sweep without an argument\n'
          'each degree has a 30-minute wall-time limit and an estimated workspace limit\n'
          'set OMP_NUM_THREADS and GRAVMAG_SPECTRAL_MAX_MIB to control resources')
    raise SystemExit(0)

# fix resource use before importing numerical libraries
os.environ.setdefault('OMP_NUM_THREADS','8')
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ.setdefault('GRAVMAG_SPECTRAL_MAX_MIB','32768')
os.environ['GRAVMAG_DIAGNOSTICS']='1'

import mpmath as mp
import numpy as np
import pandas as pd
import scipy
from gravmagpy import block,sphere_model,run_sphere_model,spherical_to_cartesian,build_fortran
import gravmagpy

root=Path(__file__).resolve().parents[1]
package=Path(gravmagpy.__file__).resolve().parents[2]
sys.path.insert(0,str(package/'tests'))
from native_support import native_library,jacobi_reference

destination=root/'diagnostics/native_high_degree'
destination.mkdir(parents=True,exist_ok=True)
try:
    physical_memory_gib=os.sysconf('SC_PHYS_PAGES')*os.sysconf('SC_PAGE_SIZE')/2**30
except (AttributeError,ValueError,OSError):
    physical_memory_gib=None


def resource_peak():
    factor=1 if sys.platform=='darwin' else 1024
    return {'process_peak_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*factor/2**20,
            'child_peak_mib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*factor/2**20}


def source_check(degree,kind):
    radius=1737.4
    altitude=60*radius/degree
    spacing=round(180/degree,6)
    latitudes=7+spacing*np.arange(-2,3)
    longitudes=-59+spacing*np.arange(-2,3)
    material={'density_kg_m3':400} if kind=='gravity' else {'magnetization_a_m':[1,-2,3]}
    source=block([6,8],[-60,-58],[2,8],**material)
    model=sphere_model([source],latitude_deg=latitudes,longitude_deg=longitudes,
                       altitude_km=altitude,radius_km=radius)
    options={'lmax':degree,'ntheta_fit':degree+1,'nphi_fit':2*degree+1,
             'refine_factor':1,'source_nlat':1,'source_nlon':1,'source_nr':1,
             'auto_mode':0,'edge_correction':0,'hybrid_mode':0,'reg_lambda':0}
    start=perf_counter()
    result=run_sphere_model(model,solver='spectral',solver_options=options,
                           output_dir=destination/f'degree_{degree}'/kind,timeout=1800)
    elapsed=perf_counter()-start
    # reconstruct the executed observation geometry from the input card precision
    cards=result['input_path'].read_text().splitlines()
    card=np.array(cards[1].split(),dtype=float)
    lat0,lon0,dlat,dlon,height,nlat,nlon=card
    latitude,longitude=np.meshgrid(lat0+dlat*np.arange(int(nlat)),lon0+dlon*np.arange(int(nlon)),indexing='ij')
    position=spherical_to_cartesian(latitude,longitude,radius+height)*1000
    center=spherical_to_cartesian(7,-59,radius-5)*1000
    volume=((radius-5)*1000)**2*np.cos(np.deg2rad(7))*np.deg2rad(2)**2*6000
    displacement=position-center
    distance=np.linalg.norm(displacement,axis=-1)
    if kind=='gravity':
        reference=-6.67430e-11*400*volume*displacement/distance[...,None]**3*1e5
    else:
        moment=volume*np.array([1.,-2.,3.])
        reference=1e-7*(3*(displacement@moment)[...,None]*displacement/distance[...,None]**5
                        -moment/distance[...,None]**3)*1e9
    error=np.linalg.norm(result['field']-reference)/np.linalg.norm(reference)
    np.savez_compressed(destination/f'degree_{degree}'/f'{kind}_reference.npz',
                        reference=reference,actual=result['field'])
    (destination/f'degree_{degree}'/kind/'stdout.txt').write_text(result['stdout'])
    if not np.isfinite(error) or error>2e-9:
        raise AssertionError(f'{kind} source error at degree {degree}: {error}')
    return {'check':f'{kind}_source','degree':degree,'relative_l2':float(error),
            'elapsed_s':elapsed,'altitude_km':float(height),'ntheta':degree+1,
            'nphi':2**int(np.ceil(np.log2(2*degree+1))),**resource_peak()}


if len(sys.argv)==2 and os.environ.get('GRAVMAG_NATIVE_WORKER')=='1':
    degree=int(sys.argv[1])
    run=destination/f'degree_{degree}'
    run.mkdir(parents=True,exist_ok=True)
    library=native_library()
    binary=build_fortran("spectral")
    records=[]
    high_precision=[]
    start_total=perf_counter()
    # include terms whose sectoral seeds underflow in an unscaled recurrence
    selections=[(degree,0,.3),(degree,1,1e-7),(degree,degree//3,np.pi/6),
                (degree,int(.43*degree),np.pi/6),(degree,degree,np.pi/2),
                (degree,max(1,degree//2),np.pi-.3)]
    for l,m,theta in selections:
        actual=np.zeros(3)
        library.native_selected_probe(l,m,float(theta),actual)
        reference=jacobi_reference(l,m,theta,80)
        # precision doubling guards the independently evaluated reference itself
        verified=jacobi_reference(l,m,theta,160)
        for index,name in enumerate(('legendre','colatitude_derivative','azimuthal_quotient')):
            scaled=abs(actual[index]-reference[index])/max(1.,abs(reference[index]))
            relative=abs(actual[index]-reference[index])/abs(reference[index]) if reference[index]!=0 else np.nan
            reference_change=abs(reference[index]-verified[index])/max(1.,abs(verified[index]))
            high_precision.append({'degree':degree,'order':m,'theta_rad':theta,'quantity':name,
                                   'actual':actual[index],'reference':reference[index],
                                   'scaled_error':scaled,'relative_error':relative,'reference_change':reference_change})
            if not np.isfinite(actual[index]) or scaled>2e-9 or reference_change>1e-25:
                raise AssertionError(f'function check failed: {high_precision[-1]}')
    pd.DataFrame(high_precision).to_csv(run/'multiprecision.csv',index=False)
    print(f'degree {degree}: independent function checks passed',flush=True)
    ntheta=degree+1
    nphi=2**int(np.ceil(np.log2(2*degree+1)))
    random=np.random.default_rng(19102026+degree)
    coefficients=random.uniform(-1.,1.,size=degree*(degree+2))
    recovered=np.empty_like(coefficients)
    statistics=np.zeros(3)
    start=perf_counter()
    library.native_roundtrip_probe(degree,ntheta,nphi,1.,coefficients,recovered,statistics)
    elapsed=perf_counter()-start
    records.append({'check':'full_bandwidth_vector_roundtrip','degree':degree,'elapsed_s':elapsed,
                    'relative_l2':statistics[2],'scaled_max_error':statistics[1],
                    'weight_sum_error':abs(statistics[0]-2),'ntheta':ntheta,'nphi':nphi,**resource_peak()})
    if not np.isfinite(statistics).all() or statistics[1]>2e-9 or statistics[2]>2e-10:
        raise AssertionError(f'roundtrip failed: {statistics}')
    del coefficients,recovered
    gc.collect()
    pd.DataFrame(records).to_csv(run/'results.csv',index=False)
    print(f'degree {degree}: full-bandwidth roundtrip passed in {elapsed:.3f} s, relative error {statistics[2]:.3e}',flush=True)
    for kind in ('gravity','magnetic'):
        record=source_check(degree,kind)
        records.append(record)
        pd.DataFrame(records).to_csv(run/'results.csv',index=False)
        print(f'degree {degree}: {kind} source passed in {record["elapsed_s"]:.3f} s, relative error {record["relative_l2"]:.3e}',flush=True)
    metadata={'degree':degree,'complete':True,'elapsed_s':perf_counter()-start_total,
              'threads':os.environ['OMP_NUM_THREADS'],'blas_threads':os.environ['OPENBLAS_NUM_THREADS'],'python':platform.python_version(),
              'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,
              'shtools_used':False,'seed':19102026+degree,
              'executable_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),**resource_peak()}
    (run/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
else:
    if len(sys.argv)>2:
        raise SystemExit('supply no argument for the resource-limited sweep or one degree for a single validation')
    native_library()
    build_fortran('spectral')
    degrees=(int(sys.argv[1]),) if len(sys.argv)==2 else (600,1200,2800,4096,5600,8192,10000,11200,14000,18000,22000)
    if any(degree<1 or degree>46339 for degree in degrees):
        raise SystemExit('degrees must be between 1 and 46339')
    attempts=[]
    existing=destination/'attempts.csv'
    if len(sys.argv)==2 and existing.exists():
        attempts=pd.read_csv(existing).to_dict('records')
        binary_hash=hashlib.sha256(build_fortran('spectral').read_bytes()).hexdigest()
        for item in attempts:
            if item['status']=='passed':
                prior=json.loads((destination/f'degree_{item["degree"]}/metadata.json').read_text())
                if prior.get('executable_sha256')!=binary_hash:
                    raise SystemExit('existing results use a different executable, rerun the full sweep')
    worker_environment={**os.environ,'GRAVMAG_NATIVE_WORKER':'1'}
    for degree in degrees:
        attempts=[item for item in attempts if item['degree']!=degree]
        ntheta=degree+1
        nphi=2**int(np.ceil(np.log2(2*degree+1)))
        estimated_bytes=60*ntheta*(degree+1)+16*degree*(degree+2)+32*(80*ntheta+80*nphi)
        if estimated_bytes>float(os.environ['GRAVMAG_SPECTRAL_MAX_MIB'])*2**20:
            attempts.append({'degree':degree,'status':'memory_limit','estimated_workspace_mib':estimated_bytes/2**20})
            break
        log=destination/f'degree_{degree}.log'
        start=perf_counter()
        with log.open('w') as stream:
            process=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),str(degree)],stdout=stream,
                                     stderr=subprocess.STDOUT,start_new_session=True,env=worker_environment)
            try:
                returncode=process.wait(timeout=1800)
                status='passed' if returncode==0 else 'failed'
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL)
                    process.wait()
                status='time_limit'
        attempts.append({'degree':degree,'status':status,'elapsed_s':perf_counter()-start,
                         'estimated_workspace_mib':estimated_bytes/2**20})
        pd.DataFrame(attempts).to_csv(destination/'attempts.csv',index=False)
        print(log.read_text(),end='',flush=True)
        print(f'degree {degree}: {status}',flush=True)
        if status!='passed':
            break
    attempts=sorted(attempts,key=lambda item:item['degree'])
    pd.DataFrame(attempts).to_csv(destination/'attempts.csv',index=False)
    complete=[item['degree'] for item in attempts if item['status']=='passed']
    tables=[pd.read_csv(destination/f'degree_{degree}'/'results.csv') for degree in complete]
    if tables:
        pd.concat(tables,ignore_index=True).to_csv(destination/'results.csv',index=False)
        pd.concat([pd.read_csv(destination/f'degree_{degree}'/'multiprecision.csv') for degree in complete],
                  ignore_index=True).to_csv(destination/'multiprecision.csv',index=False)
    (destination/'limits.json').write_text(json.dumps({'per_degree_timeout_s':1800,
        'workspace_limit_mib':float(os.environ['GRAVMAG_SPECTRAL_MAX_MIB']),'threads':os.environ['OMP_NUM_THREADS'],
        'physical_memory_gib':physical_memory_gib,'logical_cpus':os.cpu_count(),'highest_completed_degree':max(complete,default=0),
        'definition':'complete independent function checks, full-bandwidth vector roundtrip, and both public source-field checks',
        'shtools_used':False,'attempts':attempts},indent=2)+'\n')
    required=degrees if len(sys.argv)==2 else (600,1200,2800)
    if not all(degree in complete for degree in required):
        raise SystemExit('required degrees did not all pass')
