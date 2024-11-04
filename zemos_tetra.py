from src import dump, data
from src.dump import TetrahedralDistributionAnalysis, Dump
import sys

print(f"Starting ${sys.argv[3]} run")
tetrahedral_analysis = TetrahedralDistributionAnalysis(bin_height=0.01)

dump = Dump(
    sys.argv[1],
    data_file=sys.argv[2],
    outfile=sys.argv[3]
    )
    # , start=2100000
dump.ingest(frame_funcs = [tetrahedral_analysis], save_output=True, show_step=True)

tetrahedral_analysis.finalise_analysis(filename=sys.argv[3])