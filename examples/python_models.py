"""construct an equivalent volume layer without manually editing input cards"""

from pathlib import Path

import numpy as np

from gravmagpy import equivalent_source_grid, sphere_model, run_sphere_model
from gravmagpy.plotting import plot_fields

output_dir = Path(__file__).resolve().parent / "output" / "python_models"
figs_dir = Path(__file__).resolve().parent / "figs" / "python_models"
magnetization = np.zeros((1, 2, 3, 3))
magnetization[..., 0] = [1, -1, 1]
sources = equivalent_source_grid(
    latitude_edges_deg=[6, 7, 8], longitude_edges_deg=[-62, -61, -60, -59],
    depth_edges_km=[5, 10], magnetization_a_m=magnetization,
)
model = sphere_model(sources, latitude_deg=np.linspace(5, 9, 31),
                     longitude_deg=np.linspace(-63, -58, 41), altitude_km=40)
result = run_sphere_model(model, solver="point_source", solver_options={"refine_factor": 2}, output_dir=output_dir)
plot_fields(result["input_path"], result["output_path"], figs_dir / "field.png")
print(f"generated {len(sources)} source cells; peak total field {result['btot_nt'].max():.4f} nt")
print(f"wrote {output_dir}")
