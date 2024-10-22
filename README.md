Python library for reading and writing lammps data files as well as xyz files for the same system. Currently developing various analyses as part of my PhD research activities, the ultimate goal is to make this a bit more fleshed out and generalised so it can be used broadly with lammps data files rather than my individual use cases. At present the reading/writing functionalities are a little sparse but the code is written to be extensible and in theory can support the entire data file spec, though not all of it has been tested, it also currently only supports atom styles atomic and full though in principle more can be added very easily. 

Currently the dump/analysis processing is limited to my individual use cases although I am slowly trying to improve this and get a better class structure going. Will ultimately probably publish this on PyPI and see what people think but that's a fair way off at the moment.

Current Limitations:
- Assumes there is a trailing space after section titles (e.g. "Atoms # full") and after 
- Doesn't store/specify units
- Discards periodic image information (WIP)
