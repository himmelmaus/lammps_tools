from src import dump, data
from src.dump import BinnedAngleAnalysis, NaiveDiffusionAnalysis, BinnedDiffusionAnalysis, Dump, TetrahedralAngleAnalysis

# angle_analysis.finalise_analysis(filename="mean_orientation_non_pol")

# diffusion_analysis = BinnedDiffusionAnalysis(
#                                     bin_height=2.5,
#                                     dt = 1e3,
#                                     tau_range = 1e6,
#                                     tau_step = 1e3
#                                     )

spce_tetrahedral_analysis = TetrahedralAngleAnalysis(bin_height=0.22)
spce_tetrahedral_analysis.C_types = ["CC"]

# angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
print("Starting SPC/E CC")
spce_dump = Dump(
    "/home/esmith/lammps_tools/tmp_water_models/spce_prod.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
    data_file= "/home/esmith/lammps_tools/tmp_water_models/spce_data.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
    )
    # , start=2100000
spce_dump.ingest(frame_funcs = [spce_tetrahedral_analysis], save_output=True, skip=2500)

spce_tetrahedral_analysis.finalise_analysis(filename="binned_q_spce_cc_full")

spce_tetrahedral_analysis = TetrahedralAngleAnalysis(bin_height=0.22)
# spce_tetrahedral_analysis.C_types = ["CC"]

# angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
print("Starting SPC/E CT")
spce_dump = Dump(
    "/home/esmith/lammps_tools/tmp_water_models/spce_prod.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
    data_file= "/home/esmith/lammps_tools/tmp_water_models/spce_data.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
    )
    # , start=2100000
spce_dump.ingest(frame_funcs = [spce_tetrahedral_analysis], save_output=True, skip=2400)

spce_tetrahedral_analysis.finalise_analysis(filename="binned_q_spce_full")

# print("Starting TIP3p CT2")
# tip3p_tetrahedral_analysis = TetrahedralAngleAnalysis(bin_height=0.5)

# # angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
# tip3p_dump = Dump(
#     "/home/esmith/lammps_tools/tmp_water_models/tip3p_prod.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
#     data_file= "/home/esmith/lammps_tools/tmp_water_models/tip3p_data.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
#     )
#     # , start=2100000
# tip3p_dump.ingest(frame_funcs = [tip3p_tetrahedral_analysis], save_output=False, max_steps=12000000, skip=5000)

# tip3p_tetrahedral_analysis.finalise_analysis(filename="binned_q_tip3p")
# # angle_analysis.finalise_analysis(filename="mean_orientation_water_pol")

# print("Starting TIP3p CC")

# tip3p_tetrahedral_analysis = TetrahedralAngleAnalysis(bin_height=0.5)
# tip3p_tetrahedral_analysis.C_types = ["CC"]

# # angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
# tip3p_dump = Dump(
#     "/home/esmith/lammps_tools/tmp_water_models/tip3p_prod.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
#     data_file= "/home/esmith/lammps_tools/tmp_water_models/tip3p_data.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
#     )
#     # , start=2100000
# tip3p_dump.ingest(frame_funcs = [tip3p_tetrahedral_analysis], save_output=False, max_steps=12000000, skip=5000)

# tip3p_tetrahedral_analysis.finalise_analysis(filename="binned_q_tip3p_cc")