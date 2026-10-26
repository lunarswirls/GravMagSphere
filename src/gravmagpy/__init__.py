"""gravmagpy python frontend for the gravmagsphere solver and inversion toolkit"""

from .data import load_lpmag_csv
from .cards import read_input, write_source_input, run_grid_model
from .forward import predict_field
from .geometry import block_volume, parameter_names, spherical_to_cartesian
from .inversion import fit_sources, save_fit
from .equivalent import fit_equivalent_sources, equivalent_source_table, predict_equivalent_field, equivalent_field_grid
from .sphere import block, polygon, equivalent_source_grid, sphere_model, read_model_input, write_model_input, run_sphere_model
from .utils import build_fortran, run_fortran
from .ravat import (along_track_differences, regional_monopole_grid, fit_ravat_l2resl1,
                    predict_ravat_field, save_ravat_fit, load_moonbroom_wake)

__all__ = [
    "load_lpmag_csv", "predict_field", "fit_sources", "save_fit", "build_fortran", "run_fortran",
    "block_volume", "parameter_names", "spherical_to_cartesian",
    "read_input", "write_source_input", "run_grid_model",
    "block", "polygon", "equivalent_source_grid", "sphere_model", "read_model_input", "write_model_input",
    "run_sphere_model", "fit_equivalent_sources",
    "equivalent_source_table", "predict_equivalent_field", "equivalent_field_grid",
    "along_track_differences", "regional_monopole_grid", "fit_ravat_l2resl1",
    "predict_ravat_field", "save_ravat_fit",
    "load_moonbroom_wake",
]
