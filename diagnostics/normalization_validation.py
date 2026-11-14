"""Check both normalized public solvers against analytical fields at degree 600."""
import os
import sys
from pathlib import Path

if sys.argv[1:]:
    if sys.argv[1:] in (["--help"], ["-h"]):
        print("usage: python diagnostics/normalization_validation.py\n"
              "check degree-600 coefficient roundtrips and magnetic/gravity source fields for options 2/3")
        raise SystemExit(0)
    raise SystemExit("this diagnostic accepts no arguments")
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root/"src"))
sys.path.insert(0, str(root/"tests"))
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import hashlib,json
from time import perf_counter
import numpy as np
from gravmagpy import block,sphere_model,run_sphere_model,spherical_to_cartesian,build_fortran
from native_support import native_library
library=native_library()
records=[]
degree=600
rng=np.random.default_rng(19102026+degree)
coefficients=rng.uniform(-1,1,size=degree*(degree+2))
for normalization in (2,3):
    recovered=np.empty_like(coefficients);statistics=np.zeros(3)
    start=perf_counter()
    library.normalized_roundtrip_probe(degree,degree+1,2048,1.,normalization,coefficients,recovered,statistics)
    assert statistics[1]<2e-9 and statistics[2]<2e-10,statistics
    records.append({'normalization':normalization,'check':'full_bandwidth_vector_projection','degree':degree,'relative_l2':float(statistics[2]),'scaled_max':float(statistics[1]),'elapsed_s':perf_counter()-start})
    for kind,material in [('magnetic',{'magnetization_a_m':[1,-2,3]}),('gravity',{'density_kg_m3':400})]:
        source=block([6,8],[-60,-58],[2,8],**material)
        altitude=60*1737.4/degree
        model=sphere_model([source],latitude_deg=7+.3*np.arange(-2,3),longitude_deg=-59+.3*np.arange(-2,3),altitude_km=altitude)
        options={'normalization':normalization,'lmax':degree,'ntheta_fit':degree+1,'nphi_fit':2*degree+1,'refine_factor':1,'source_nlat':1,'source_nlon':1,'source_nr':1,'auto_mode':0,'edge_correction':0,'hybrid_mode':0,'reg_lambda':0}
        start=perf_counter();result=run_sphere_model(model,solver='gauss_legendre',solver_options=options)
        lat,lon=np.meshgrid(result['latitude_deg'],result['longitude_deg'],indexing='ij')
        observation=spherical_to_cartesian(lat,lon,1737.4+altitude)*1000
        center=spherical_to_cartesian(7,-59,1732.4)*1000
        volume=(1732.4*1000)**2*np.cos(np.deg2rad(7))*np.deg2rad(2)**2*6000
        delta=observation-center;distance=np.linalg.norm(delta,axis=-1)
        if kind=='magnetic':
            moment=volume*np.array([1,-2,3]);reference=100*(3*(delta@moment)[...,None]*delta/distance[...,None]**5-moment/distance[...,None]**3)
        else:
            reference=-6.67430e-6*400*volume*delta/distance[...,None]**3
        relative=np.linalg.norm(result['field']-reference)/np.linalg.norm(reference)
        assert relative<2e-9,relative
        records.append({'normalization':normalization,'check':kind+'_source','degree':degree,'relative_l2':float(relative),'elapsed_s':perf_counter()-start})
report={'threads':int(os.environ['OMP_NUM_THREADS']),'degree':degree,'binary_sha256':hashlib.sha256(build_fortran('gauss_legendre').read_bytes()).hexdigest(),'checks':records}
(root/'diagnostics/normalization_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
