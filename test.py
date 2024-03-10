from src import dump, data
from src.dump import BinnedAngleAnalysis, NaiveDiffusionAnalysis, BinnedDiffusionAnalysis, Dump

# analysis = BinnedAngleAnalysis(bin_height=0.25, wrap_angle=True)
analysis = BinnedDiffusionAnalysis(
                                    bin_height=1,
                                    dt = 2e3,
                                    tau_range = 2.5e6,
                                    tau_step = 2e3
                                    )
dump = Dump(
    "/home/esmith/water_gold_slab/polarisable/dpg_poster/111_water_large/trajectories/water_gold_nvt_unwrapped.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
    data_file= "/home/esmith/water_gold_slab/polarisable/dpg_poster/111_water_large/data_files/gold_water_npt.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
    )
dump.ingest(frame_funcs = [analysis], save_output=False, skip=2e3, start=3e6)

analysis.finalise_analysis()