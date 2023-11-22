from typing import Any
import numpy as np
from .section import Section
from .helpers import num_str, snake_case
import re


class Data:
    hkeywords = [
        ["atoms", "atoms"],
        ["ellipsoids", "ellipsoids"],
        ["lines", "lines"],
        ["triangles", "triangles"],
        ["bodeies", "bodies"],
        ["bonds", "bonds"],
        ["angles", "angles"],
        ["dihedrals", "dihedrals"],
        ["impropers", "impropers"],
        ["atom_types", "atom types"],
        ["bond_types", "bond types"],
        ["angle_types", "angle types"],
        ["dihedral_types", "dihedral types"],
        ["improper_types", "improper types"],
        ["x_bounds", "xlo xhi"],
        ["y_bounds", "ylo yhi"],
        ["z_bounds", "zlo zhi"],
        ["tilts", "xy xz yz"],
    ]

    skeywords = [
        ["Masses", "atom_types"],
        ["Atoms", "atoms"],
        ["Ellipsoids", "ellipsoids"],
        ["Lines", "lines"],
        ["Triangles", "triangles"],
        ["Bodies", "bodies"],
        ["Bonds", "bonds"],
        ["Angles", "angles"],
        ["Dihedrals", "dihedrals"],
        ["Impropers", "impropers"],
        ["Velocities", "atoms"],
        ["Pair Coeffs", "atom_types"],
        ["Bond Coeffs", "bond_types"],
        ["Angle Coeffs", "angle_types"],
        ["Dihedral Coeffs", "dihedral_types"],
        ["Improper Coeffs", "improper_types"],
        ["BondBond Coeffs", "angle_types"],
        ["BondAngle Coeffs", "angle_types"],
        ["MiddleBondTorsion Coeffs", "dihedral_types"],
        ["EndBondTorsion Coeffs", "dihedral_types"],
        ["AngleTorsion Coeffs", "dihedral_types"],
        ["AngleAngleTorsion Coeffs", "dihedral_types"],
        ["BondBond13 Coeffs", "dihedral_types"],
        ["AngleAngle Coeffs", "improper_types"],
        ["Molecules", "atoms"],
    ]

    def __init__(self):
        self.hvars = list(np.array(self.hkeywords).T[0])
        self.hnames = list(np.array(self.hkeywords).T[1])

        self.snames = list(np.array(self.skeywords).T[0])
        self.slength = list(np.array(self.skeywords).T[1])
        self.sections = []
        self.atom_type_labels = None
        self.units = None
        
    def __getattr__(self, item):
        """
        Implementing custom setter/getter logic to allow for sections to be accessible by name,
        while still stored in an array. I'm not sure this is the best way to do it but will have
        a ponder. These need to be prepended with section_ so as to avoid confusion between e.g.
        the Atoms section and the number of atoms header value.
        """
        try:
            section_id = re.match(r"section_(.+)", item).group(1)
            return self.sections[self.section_ids.index(section_id)]
        except AttributeError:
            return super(Data, self).__getattribute__(item)
        
    def __setattr__(self, item, value):
        section_id = re.match(r"section_(.+)", item)
        if section_id:
            section_id = section_id.group(1)
            self.sections[self.section_ids.index(section_id)] = value
        else:
            return super(Data, self).__setattr__(item, value)
        
        
    @property
    def section_ids(self):
        return list(map(lambda x: x.section_id, self.sections))

    @classmethod
    def read_file(cls, filename):
        instance = cls()
        with open(filename, "r") as infile:
            lines = infile.readlines()
            instance.process_file(lines)
            return instance

    def process_file(self, lines):
        for i, line in enumerate(lines):
            if any(header_name in line for header_name in self.hnames):
                self.process_header(line.strip().split(" "))
                continue
            if re.search(r"units = (.+)\n", line):
                self.units = re.search(r"units = (.+)\n", line).group(1)
            if any(section_name in line for section_name in self.snames) or "Atom Type Labels" in line:
                break

        while i < len(lines):
            line = lines[i]
            if any(section_name in line for section_name in self.snames):
                title_style = line.strip("\n").split(" # ")
                section_lines = lines[
                    i + 2 : i + int(getattr(self,self.slength[self.snames.index(title_style[0].strip())],)) + 2
                ]
                self.sections.append(
                    Section(
                        section_lines,
                        title_style[0],
                        style=None if len(title_style) == 1 else title_style[1],
                    )
                )
                i = (
                    i + int(getattr(self,self.slength[self.snames.index(title_style[0].strip())],)) + 2
                )

            # requires custom logic as type labels are non numeric
            if "Atom Type Labels" in line:
                section_lines = lines[i + 2 : i + int(getattr(self, "atom_types")) + 2]
                _, self.atom_type_labels = np.array(
                    list(map(lambda x: x.strip("\n").split(" "), section_lines))
                ).T
                
                i = i + int(getattr(self, "atom_types")) + 2

            else:
                i = i + 1

        return self

    def process_header(self, tokens: list[str]):
        # Most header values, with a single value and a label with a space in it
        if len(tokens) <= 3:
            setattr(
                self,
                self.hvars[self.hnames.index(" ".join(tokens[1:]))],
                float(tokens[0]),
            )

        # coord lo/hi values
        elif len(tokens) == 4:
            setattr(
                self,
                self.hvars[self.hnames.index(" ".join(tokens[2:]))],
                np.array([float(tokens[0]), float(tokens[1])]),
            )

        # tilt factors
        elif len(tokens) == 6:
            setattr(
                self,
                self.hvars[self.hnames.index(" ".join(tokens[3:]))],
                np.array(list(map(float, tokens[:3]))),
            )
            
    def add_section(self, lines, title, style=None):
        if title in self.snames and snake_case(title) not in self.section_ids:
            self.sections.append(
                Section(
                    lines,
                    title,
                    style=style
                )
            )
            return self.sections[-1]
        else:
            if title in self.snames:
                raise NameError(f"Section {title} already exists")
            raise NameError(f"{title} is not a valid section title")

    def write_file(self, file="out.data", sort=True):
        outlines = [
            f"LAMMPS Data File written via lammps_tools, Elspeth Smith, Ruhr Universitaet Bochum{f', units = {self.units}' if self.units else ''}\n\n"
        ]
        for i, field in enumerate(self.hvars):
            if getattr(self, field, None) is not None:
                outlines.append(f"{num_str(getattr(self, field))} {self.hnames[i]}\n")
        outlines.append("\n\n")
        
        if getattr(self, "atom_type_labels", None) is not None:
            lines = list(map(lambda x: f"{x[0]+1} {x[1]}\n", enumerate(self.atom_type_labels)))
            outlines.append("Atom Type Labels\n\n")
            outlines.extend(lines)
            outlines.append("\n")

        for section in self.sections:
            outlines.extend(section.print_lines(sort=sort))
            outlines.append("\n")

        with open(file, "w") as outfile:
            outfile.writelines(outlines)
        return

    def write_xyz(self, file="out.xyz"):
        atoms = [section for section in self.sections if section.title=="Atoms"][0]
        outlines = [f"{len(atoms.lines)}\n"]
        outlines.extend(
            [
                f"{self.atom_type_labels[atom.type_id-1]} {' '.join(map(str, atom.coords))}\n"
                for atom in atoms.lines
            ]
        )

        with open(file, "w") as outfile:
            outfile.writelines(outlines)
        return
