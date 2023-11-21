from data import Data

data = Data.read_file("/home/esmith/lammps_tools/gold_test.data")


data.write_file()
data.write_xyz()