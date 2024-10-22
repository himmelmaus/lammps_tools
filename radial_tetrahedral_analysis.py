from src import dump, data
from src.dump import BinnedAngleAnalysis, NaiveDiffusionAnalysis, BinnedDiffusionAnalysis, Dump, TetrahedralAngleAnalysis
import threading


def get_tetrahedral_analysis(
                        analysis_name,
                        dump_file,
                        data_file = None,
                        C_types = None,
                        O_types = None,
                        bin_height=0.5,
                        ingest_kwargs = {
                            "save_output": False,
                            "skip": 2500
                        },
                        outfile_name = "binned_q"
                        ):
    def return_func():
        print(f"Starting {analysis_name}")
        tetrahedral_analysis = TetrahedralAngleAnalysis(bin_height=bin_height)

        if C_types is not None:
            tetrahedral_analysis.C_types = C_types
        
        if O_types is not None:
            tetrahedral_analysis.O_types = O_types


        # angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
        dump = Dump(
            dump_file,
            data_file = data_file
            )
            # , start=2100000
        dump.ingest(frame_funcs = [tetrahedral_analysis], **ingest_kwargs)

        tetrahedral_analysis.finalise_analysis(filename=outfile_name)
        print(f"Finished {analysis_name}")

        return True
    return return_func

tip3p_ct2 = get_tetrahedral_analysis(
    "TIP3P CT2",
    "/home/esmith/lammps_tools/tmp_water_models/tip3p_prod.xyz",
    "/home/esmith/lammps_tools/tmp_water_models/tip3p_data.data",
    C_types = ["CT2"],
    bin_height=0.25,
    outfile_name = "q_outfiles/tip3p_binned_q_ct2"
)

tip3p_cc = get_tetrahedral_analysis(
    "TIP3P CC",
    "/home/esmith/lammps_tools/tmp_water_models/tip3p_prod.xyz",
    "/home/esmith/lammps_tools/tmp_water_models/tip3p_data.data",
    C_types = ["CC"],
    bin_height=0.25,
    outfile_name = "q_outfiles/tip3p_binned_q_cc"
)

spce_ct2 = get_tetrahedral_analysis(
    "SPC/E CT2",
    "/home/esmith/lammps_tools/tmp_water_models/spce_prod_full.xyz",
    "/home/esmith/lammps_tools/tmp_water_models/spce_data.data",
    C_types = ["CT2"],
    bin_height=0.25,
    outfile_name = "q_outfiles/spce_binned_q_ct2"
)

spce_cc = get_tetrahedral_analysis(
    "SPC/E CC",
    "/home/esmith/lammps_tools/tmp_water_models/spce_prod_full.xyz",
    "/home/esmith/lammps_tools/tmp_water_models/spce_data.data",
    C_types = ["CC"],
    bin_height=0.25,
    outfile_name = "q_outfiles/spce_binned_q_cc"
)

tip3p_ct2_thread = threading.Thread(target=tip3p_ct2)
tip3p_cc_thread = threading.Thread(target=tip3p_cc)

spce_ct2_thread = threading.Thread(target=spce_ct2)
spce_cc_thread = threading.Thread(target=spce_cc)

tip3p_ct2_thread.start()
tip3p_cc_thread.start()

spce_ct2_thread.start()
spce_cc_thread.start()




# print("Starting TIP3p CC")

# tip3p_tetrahedral_analysis = TetrahedralAngleAnalysis(bin_height=0.5)
# tip3p_tetrahedral_analysis.C_types = ["CC"]

# # angle_analysis = BinnedAngleAnalysis(bin_height=2.5)
# tip3p_dump = Dump(
#     "/home/esmith/lammps_tools/tmp_water_models/tip3p_prod.xyz", #"/home/esmith/water_gold_slab/polarisable/111_water/trajectories/water_gold_nvt_unwrapped.xyz", 
#     data_file= "/home/esmith/lammps_tools/tmp_water_models/tip3p_data.data" # "/home/esmith/water_gold_slab/polarisable/111_water/data_files/gold_water_npt.data"
#     )
#     # , start=2100000
# tip3p_dump.ingest(frame_funcs = [tip3p_tetrahedral_analysis], save_output=False, max_steps=15000000, skip=2500)

# tip3p_tetrahedral_analysis.finalise_analysis(filename="binned_q_tip3p_cc")