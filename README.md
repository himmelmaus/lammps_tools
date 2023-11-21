Python library for reading and writing lammps data files, currently no functionality beyond this. More advanced utilities/refactoring to come.

Current Limitations:
- Assumes there is a trailing space after section titles (e.g. "Atoms # full") and after 
- Doesn't store/specify units
- Discards periodic image information
- ~~TODO: make sections directly accessible by name~~