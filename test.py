from src import dump, data
from src.dump import BinnedAngleAnalysis, NaiveDiffusionAnalysis, BinnedDiffusionAnalysis, Dump

# angle_analysis.finalise_analysis(filename="mean_orientation_non_pol")

diffusion_analysis = BinnedDiffusionAnalysis(
                                    bin_height=2.5,
                                    dt = 1e3,
                                    tau_range = 1e6,
                                    tau_step = 1e3
                                    )

# angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
dump = Dump(
    "/home/esmith/water_gold_slab/polarisable/dpg_poster/111_water_large/trajectories/water_gold_nvt_unwrapped.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
    data_file= "/home/esmith/water_gold_slab/polarisable/dpg_poster/111_water_large/data_files/gold_water_nvt_10_ns.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
    )
dump.ingest(frame_funcs = [diffusion_analysis], save_output=False, skip=1e3, start=5e6)

diffusion_analysis.finalise_analysis(filename="binned_pol_rabe")
# angle_analysis.finalise_analysis(filename="mean_orientation_water_pol")