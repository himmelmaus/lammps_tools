from typing import Any
from itertools import product
import multiprocessing as mp
from datetime import datetime

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from src.dump import Dump
from src.atom import Atom


class Analysis:

    """
    Example analysis script to calculate the RMS angle of water molecules relative to the x-axis
    """
    
    unit_z = [0,0,1]
    
    def __init__(self, bin_height=0.5, fixed_bins=True, **kwargs):
        self.bin_height = bin_height
        self.initial_bin_height = bin_height
        self.mol_atoms = None
        self.fixed_bins = fixed_bins
        self.bins = None

    
    def __call__(self, old_atoms, new_atoms, *, dump):

        if self.bins is None or not self.fixed_bins:
            self.nbins = int(np.ceil(np.abs(dump.zbounds[1] - dump.zbounds[0])/self.bin_height))
            self.bin_height = np.abs(dump.zbounds[1] - dump.zbounds[0])/self.nbins
            self.bins = np.linspace(dump.zbounds[0], dump.zbounds[1], self.bin_height)

        if self.mol_atoms is None:
            self.set_mol_atoms(dump)

        angles = []
        for atoms in self.mol_atoms.values():

            angles.append(self.get_mol_angle(atoms))

        print(dump.current_step, dump.n_atoms, np.mean(angles), max(angles))
        return np.mean(angles)
    
    def set_mol_atoms(self, dump):
        self.mol_atoms = {}
        O_H_atoms = [a.atom_id for a in dump.atoms if a.type_id in ["O", "H"]]
        for mol_id, atoms in dump.molecules.items():
            # Water specific sorting so that the oxygen is always in the first position
            # this will also break if the analysis is applied to multiple dumps
            if any(a not in O_H_atoms for a in atoms) or mol_id == 0:
                continue
                

            self.mol_atoms[mol_id] = sorted(filter(lambda a: a.atom_id in atoms, dump.atoms), key=lambda a: 0 if a.type_id == "O" else 1)

    def get_mol_vector(self, atoms):
        H_midpoint = (atoms[2].coords + atoms[1].coords)/2
        orientation_vector = np.copy(H_midpoint - atoms[0].coords)
        orientation_vector /= np.linalg.norm(orientation_vector) # normalising

        return orientation_vector
    
    def get_mol_angle(self, atoms):

        orientation_vector = self.get_mol_vector(atoms)

        # Radians!!!!!
        angle = np.arccos(np.dot(orientation_vector, self.unit_z))

        if self.wrap_angle:
            angle = np.pi - angle if angle > np.pi/2 else angle # CAST

        return angle
    
    def get_bin_index(self, coords, zbounds):

        if coords[2] >= zbounds[1]:
            bin_index = int(np.floor((coords[2]%zbounds[1] + zbounds[0])/self.bin_height))
        elif coords[2] < zbounds[0]:
            box_height = zbounds[1] - zbounds[0]
            wrapped_coord = coords[2] + (1 + np.floor(np.abs(coords[2])/box_height))*box_height
            bin_index = int(np.floor(wrapped_coord/self.bin_height))
        else:
            bin_index = int(np.floor(coords[2]/self.bin_height))

        return bin_index
    
class BinnedAngleAnalysis(Analysis):

    def __init__(self, bin_height=0.5, fixed_bins=True, wrap_angle=False):

        super(BinnedAngleAnalysis, self).__init__(bin_height=bin_height, fixed_bins=fixed_bins, wrap_angle=wrap_angle)
        self.wrap_angle = wrap_angle

        if not fixed_bins:
            raise NotImplementedError("Support for variable bin height has not yet been implemented for this analysis.")

        # Rather than calculating a rolling average or similar, keep track of cumulative sum and number of entries
        self.cum_bin_angles = None 
        self.n_bin_entries = None


    def __call__(self, old_atoms, new_atoms, *, dump):

        if self.bins is None:
            
            # Shrink the bin size slightly to accommodate equally sized bins
            self.nbins = int(np.ceil(np.abs(dump.zbounds[1] - dump.zbounds[0])/self.bin_height)) + 1
            self.bin_height = np.abs(dump.zbounds[1] - dump.zbounds[0])/(self.nbins-1)

            self.bins = np.linspace(dump.zbounds[0], dump.zbounds[1], self.nbins)

            self.cum_bin_angles = np.zeros((self.nbins,))
            self.n_bin_entries = np.zeros((self.nbins,))

        if self.mol_atoms is None:
            self.set_mol_atoms(dump)

        
        for atoms in self.mol_atoms.values():

            if atoms[0].coords[2] >= dump.zbounds[1]:
                bin_index = int(np.floor((atoms[0].coords[2]%dump.zbounds[1] + dump.zbounds[0])/self.bin_height))
            elif atoms[0].coords[2] < dump.zbounds[0]:
                box_height = dump.zbounds[1] - dump.zbounds[0]
                wrapped_coord = atoms[0].coords[2] + (1 + np.floor(np.abs(atoms[0].coords[2])/box_height))*box_height
                bin_index = int(np.floor(wrapped_coord/self.bin_height))
            else:
                bin_index = int(np.floor(atoms[0].coords[2]/self.bin_height))
                
            self.cum_bin_angles[bin_index] += self.get_mol_angle(atoms)
            self.n_bin_entries[bin_index] += 1
        
        return
    
    def finalise_analysis(self, save=True, filename="mean_orientation_binned"):

        mean_orientations = self.cum_bin_angles/self.n_bin_entries
        mean_orientations = np.nan_to_num(mean_orientations, posinf=0, neginf=0)
        mean_orientations_deg = np.degrees(mean_orientations)

        plt.plot(self.bins, mean_orientations_deg)
        plt.xlabel("Z Coordinate(Å)")
        if self.wrap_angle:
            plt.ylabel("Mean lowest angle between molecules and z axis (Degrees)")
        else:
            plt.ylabel("Mean angle between molecules and positive z axis (Degrees)")
        
        np.savetxt(f"{filename}.txt", (self.bins, mean_orientations_deg))
        plt.savefig(f"{filename}.png")
        plt.cla()
        plt.clf()
        return
    
class CitrateAngleAnalysis(Analysis):
    
    """
    For finding the angle made between citrate molecules and the negative z axis
    """
    
    unit_z = [0,0,1]
    
    def __init__(self, bin_height=0.5, fixed_bins=True, wrap_angle=False):

        super(CitrateAngleAnalysis, self).__init__(bin_height=bin_height, fixed_bins=fixed_bins, wrap_angle=wrap_angle)
        self.wrap_angle = wrap_angle

        if not fixed_bins:
            raise NotImplementedError("Support for variable bin height has not yet been implemented for this analysis.")

        # Rather than calculating a rolling average or similar, keep track of cumulative sum and number of entries
        self.cum_bin_angles = None 
        self.n_bin_entries = None



    
    def __call__(self, old_atoms, new_atoms, *, dump):

        if self.bins is None or not self.fixed_bins:
            self.nbins = int(np.ceil(np.abs(dump.zbounds[1] - dump.zbounds[0])/self.bin_height))
            self.bin_height = np.abs(dump.zbounds[1] - dump.zbounds[0])/self.nbins
            self.bins = np.linspace(dump.zbounds[0], dump.zbounds[1], self.bin_height)

        if self.mol_atoms is None:
            self.set_mol_atoms(dump)

        angles = []
        for atoms in self.mol_atoms.values():

            angles.append(self.get_mol_angle(atoms))

        print(dump.current_step, dump.n_atoms, np.mean(angles), max(angles))
        return np.mean(angles)
    
    def set_mol_atoms(self, dump):
        self.mol_atoms = {}
        CC_atoms = [a.atom_id for a in dump.atoms if a.type_id == "CC"]
        for mol_id, atoms in dump.molecules.items():
            # Water specific sorting so that the oxygen is always in the first position
            # this will also break if the analysis is applied to multiple dumps
            cc_atom_indices = np.where([atom.atom_id in CC_atoms for atom in atoms])
            # mol_cc_atoms = atoms[cc_atom_indices]
                

            self.mol_atoms[mol_id] = [atoms[np.max(cc_atom_indices)], atoms[np.min(cc_atom_indices)]]# [sorted(filter(lambda a: a.atom_id in atoms, dump.atoms), key=lambda a: a.atom_id)]

    def get_mol_vector(self, atoms):
        # H_midpoint = (atoms[2].coords + atoms[1].coords)/2
        orientation_vector = np.copy(atoms[1].coords - atoms[0].coords)
        orientation_vector /= np.linalg.norm(orientation_vector) # normalising

        return orientation_vector
    
    def get_mol_angle(self, atoms):

        orientation_vector = self.get_mol_vector(atoms)

        # Radians!!!!!
        angle = np.arccos(np.dot(orientation_vector, self.unit_z))

        if self.wrap_angle:
            angle = np.pi - angle if angle > np.pi/2 else angle # CAST

        return angle
    
    def get_bin_index(self, coords, zbounds):

        if coords[2] >= zbounds[1]:
            bin_index = int(np.floor((coords[2]%zbounds[1] + zbounds[0])/self.bin_height))
        elif coords[2] < zbounds[0]:
            box_height = zbounds[1] - zbounds[0]
            wrapped_coord = coords[2] + (1 + np.floor(np.abs(coords[2])/box_height))*box_height
            bin_index = int(np.floor(wrapped_coord/self.bin_height))
        else:
            bin_index = int(np.floor(coords[2]/self.bin_height))

        return bin_index

class NaiveDiffusionAnalysis(Analysis):
    
    oxygens = None
    oxygen_indices = None
    original_positions = None
    
    
    def __call__(self, old_atoms, new_atoms, *, dump: Dump):
        
        # print(dump.current_step)
        self.dump = dump
        
        if self.oxygen_indices is None:
            self.oxygens = np.fromiter(filter(lambda a: a.type_id == "O", dump.atoms), dtype=Atom)
            self.oxygen_indices = np.fromiter(map(lambda a: a.atom_id, self.oxygens), dtype=int)
            self.original_coords = np.fromiter(map(lambda a: np.array(a.coords)[:2], self.oxygens), dtype=(np.float64,2))

        if self.bins is None:
            
            # Shrink the bin size slightly to accommodate equally sized bins
            self.nbins = int(np.ceil(np.abs(dump.zbounds[1] - dump.zbounds[0])/self.bin_height))
            self.bin_height = np.abs(dump.zbounds[1] - dump.zbounds[0])/self.nbins

            self.bins = np.linspace(dump.zbounds[0], dump.zbounds[1], self.nbins)

            self.bin_msds = np.zeros((self.nbins, 0))
            self.timesteps = []


        bin_displacements = np.zeros(self.nbins)
        atom_count = np.zeros(self.nbins)
        self.timesteps.append(dump.current_step)

        for i, atom in enumerate(self.oxygens):

            # 2D displacement in the x-y plane
            disp = atom.coords[:2] - self.original_coords[i]
            square_disp = np.linalg.norm(disp)**2

            bin_index = self.get_bin_index(atom.coords, dump.zbounds)

            bin_displacements[bin_index] += square_disp
            atom_count[bin_index] += 1

        msds = (bin_displacements/atom_count)
        self.bin_msds = np.append(self.bin_msds, np.array([msds]).T, axis=1)
        
        return
    
    def finalise_analysis(self, save=True, filename="binned_msd_plot.png"):

        plt.plot(self.timesteps, self.bin_msds.T, label=list(map(int, self.bins)))
        plt.legend()
        plt.tight_layout()
        plt.savefig(filename)
        plt.cla()
        plt.clf()
        
        ms = np.array([])
        cs = np.array([])
        
        straight_line = lambda x, m, c: m*x + c
        
        fit_timesteps = np.array(self.timesteps)*self.dump.timestep
        
        for msd in self.bin_msds:
            
            nan_filter = (~(np.isnan(msd)) & ~(np.isinf(msd)) & ~(np.isclose(msd, 0)))
            
            if not nan_filter.any():
                ms = np.append(ms, 0)
                cs = np.append(cs, 0)
                continue
            
            timesteps = np.array(fit_timesteps)[nan_filter]
            msd = msd[nan_filter]
            
            
            popt, pcov = curve_fit(straight_line, timesteps, msd)
            
            print(popt, pcov)
            
            ms = np.append(ms, popt[0])
            cs = np.append(cs, popt[1])
        
        
        plot_timesteps = np.repeat([[min(fit_timesteps), max(fit_timesteps)]], self.nbins, axis=0)
        fit_lines = straight_line(plot_timesteps, ms.reshape((self.nbins, 1)), cs.reshape((self.nbins, 1)))

        plt.plot(plot_timesteps[0], fit_lines.T)
        plt.tight_layout()
        plt.savefig("msd_test.png")
        plt.cla()
        plt.clf()
        
        # convert to diffusion coefficients in units of (m^2)/s
        diffusion_coefficients = (ms/4)*(10**-5)
        
        plt.plot(self.bins, diffusion_coefficients)
        plt.xlabel("Bin Z Coordinate (Å)")
        plt.ylabel("Coefficient of Diffusion $m^{2}s^{-1}$")
        plt.tight_layout()
        plt.savefig("binned_diffusion_coefficient.png")
        
        return
    
class BinnedDiffusionAnalysis(Analysis):
    
    oxygens = None
    oxygen_indices = None
    n_oxygens = None
    nbins = None
    bin_msds = None
    timesteps = None
    n_threads = 4
    
    buffer_coords = None
    buffer_frames = None
    buffer_indices = None

    dt = None
    tau_range = None
    tau_step = None
    tau_bin_pairs = None
    
    ntau = None
    cum_bin_msd = None 
    cum_bin_entries = None 
    
    def __init__(self,
                 bin_height=0.5,
                 fixed_bins=True,
                 dt = 10000,
                 tau_range = 100000,
                 tau_step = 1000,
                 **kwargs):
        
        super().__init__(bin_height, fixed_bins, **kwargs)
        
        self.dt=dt
        self.tau_range = tau_range
        self.tau_step = tau_step
        self.ntau = int(self.tau_range/self.tau_step)
        
    
    def __call__(self, old_atoms, new_atoms, *, dump: Dump):
        
        """
        Brief explainer on the methodology:
        
        Using the method for binned diffusion coefficient calculation defined in Ntim & Sulpizi, 2020. This uses a rolling window from which to calculate
        MSDs with different steps for window width and gradations in the window start time. I achieve this by recording oxygen coordinates at each frame
        which corresponds to a rolling window start, as well as a mask array which records which bin each oxygen is at the specified timestep. When actually
        calculating the MSD from this, I combine the oxygen masks from each respective step to get the list of coordinates for the relevant calculation
        """

        if self.bins is None:
            
            # Shrink the bin size slightly to accommodate equally sized bins
            self.nbins = int(np.ceil(np.abs(dump.zbounds[1] - dump.zbounds[0])/self.bin_height))+1
            self.bin_height = np.abs(dump.zbounds[1] - dump.zbounds[0])/(self.nbins-1)

            self.bins = np.linspace(dump.zbounds[0], dump.zbounds[1], self.nbins)

            self.timesteps = []
            
            if dump.skip:
                assert self.dt % dump.skip == 0
            else:
                assert self.dt % dump.timestep == 0
                
            
        
        if self.oxygen_indices is None:
            self.oxygens = np.fromiter(filter(lambda a: a.type_id == "O", dump.atoms), dtype=Atom)
            self.oxygen_indices = np.fromiter(map(lambda a: a.atom_id, self.oxygens), dtype=int)
            self.n_oxygens = len(self.oxygen_indices)

            # frames will be buffered so as to minimise RAM usage, god I fucking love vectorisation
            self.buffer_frames = np.array([])
            self.buffer_indices = np.ndarray((0, self.nbins, self.n_oxygens), dtype=np.int8)
            self.buffer_coords = np.ndarray((0, self.n_oxygens, 3))
            
            self.cum_bin_msd = np.zeros((self.nbins, self.ntau))
            self.cum_bin_entries = np.zeros((self.nbins, self.ntau), dtype=np.int64)
            
            # pre computing this and saving
            self.tau_bin_pairs = np.array(list(product(range(self.ntau), range(self.nbins))))
            
        
        # not doing any actual processing yet, just gathering information in the correct
            
        # need to reshape the array in this cursed way to do np.append properly
        self.buffer_frames = np.append(self.buffer_frames, dump.current_step)
        oxygen_coords = np.reshape(
                                    np.fromiter(
                                        map(lambda a: np.array(a.coords), self.oxygens), dtype=(np.float64,3)
                                    ), 
                                    (1, self.n_oxygens, 3)
                                )
        self.buffer_coords = np.append(self.buffer_coords, oxygen_coords, axis=0)

        # given we don't need to really worry about periodic coords in the z dimension
        # given the gold slabs preventing periodicity, I think it should be possible
        # to get a list of gold atom indices by bin using just vector operations
        # but I need to go to rewe right now... 

        # okay never mind I do need to worry about periodic coords
        
        
        
        # if coords[2] >= zbounds[1]:
        #     bin_index = int(np.floor((coords[2]%zbounds[1] + zbounds[0])/self.bin_height))
        # elif coords[2] < zbounds[0]:
        #     box_height = zbounds[1] - zbounds[0]
        #     wrapped_coord = coords[2] + (1 + np.floor(np.abs(coords[2])/box_height))*box_height
        #     bin_index = int(np.floor(wrapped_coord/self.bin_height))
        
        zlo = min(dump.zbounds)
        zhi = max(dump.zbounds)
        zs = np.copy(oxygen_coords[0,:, 2])

        zs[zs > zhi] = zs[zs > zhi]%zhi + zlo
        zs[zs < zlo] = zs[zs < zlo] + (1+np.floor(np.abs(zs[zs<zlo])/(zhi-zlo)))*(zhi-zlo)

        oxygen_bin_indices = np.floor(zs/self.bin_height).astype(int)

    
        # oxygen_bin_indices = np.zeros(zs.shape)
        
        # oxygen_bin_indices[(zs >= dump.zbounds[0]) & (zs < dump.zbounds[1])] = np.floor(zs[(zs >= dump.zbounds[0]) & (zs < dump.zbounds[1])]/self.bin_height)
        
        # box_height = dump.zbounds[1] - dump.zbounds[0]
        # # removed a 1+ here
        # oxygen_bin_indices[(oxygen_bin_indices < dump.zbounds[0])] = (zs[(oxygen_bin_indices < dump.zbounds[0])] + (1 + np.floor(np.abs(zs[(oxygen_bin_indices < dump.zbounds[0])])/box_height))*box_height/self.bin_height).astype(int)
        
        # oxygen_bin_indices[(oxygen_bin_indices >= dump.zbounds[0])] = np.floor((zs[(oxygen_bin_indices >= dump.zbounds[0])]%dump.zbounds[1] + dump.zbounds[0])/self.bin_height)
        
        
        # oxygen_bin_indices = oxygen_bin_indices.astype(int)

        # create a mask of self.nbins x self.n_oxygens which will hopefully allow for easy 
        # combination of values between steps and thus easy retrieval of coords
        
        # might need to change this to just a list of which indices are in which bin at a given time
        # but then will run into issues of an inhomo array and ugh
        coords_mask = np.zeros((self.nbins, self.n_oxygens), dtype=np.int8)
        coords_mask[oxygen_bin_indices, np.arange(self.n_oxygens)] = 1

        self.buffer_indices = np.append(self.buffer_indices, [coords_mask], axis=0)
        
        
        # onto actual data processing now
        # this checks whether the time difference between the latest gathered and earliest gathered frame has reached a full tau window
        # if it has, we will discard the earliest one after calculating the msd windows starting from it
        if len(self.buffer_frames) > self.tau_range/self.dt:
            
            # with mp.Pool(self.n_threads) as pool:
            #     for ret in pool.imap_unordered(self.bin_msd, self.tau_bin_pairs):
            #         print(ret[:1])
            #         self.cum_bin_msd[ret[0], ret[1]] += ret[2]
            #         self.cum_bin_entries[ret[0], ret[1]] += ret[3]
            
            for i in range(1,self.ntau):  
                
                oxygens_filter = (self.buffer_indices[0]) & (self.buffer_indices[i])
                
                relative_vectors = self.buffer_coords[i,:,:2] - self.buffer_coords[0,:,:2]
                
                squared_disp = np.linalg.norm(relative_vectors, axis=1)**2
                
                
                self.cum_bin_msd[:,i] += np.sum(squared_disp*oxygens_filter, axis=1)
                self.cum_bin_entries[:,i] += np.sum(oxygens_filter, axis=1)
                
                # pretty sure I could vectorise this to remove the inner loop
                # for j in range(self.nbins):
                #     oxygens_filter = (self.buffer_indices[0][j]) & (self.buffer_indices[i][j])
                    
                #     # discard z coord
                #     relative_vectors = (self.buffer_coords[i][oxygens_filter][:,:2] - self.buffer_coords[0][oxygens_filter][:,:2])
                    
                #     squared_disp = np.linalg.norm(relative_vectors, axis=1)**2
                    
                #     self.cum_bin_msd[j][i] += np.sum(squared_disp)
                #     self.cum_bin_entries[j][i] += len(squared_disp)
            
            # remove the earliest frame as we're now done with that as a reference
            self.buffer_frames = np.delete(self.buffer_frames, 0, axis=0)
            self.buffer_indices = np.delete(self.buffer_indices, 0, axis=0)
            self.buffer_coords = np.delete(self.buffer_coords, 0, axis=0)
            
    def bin_msd(self, tau_bin):
        
        tau_index = tau_bin[0]
        bin_index = tau_bin[1]
        
        oxygens_filter = (self.buffer_indices[0][bin_index]) & (self.buffer_indices[tau_index][bin_index])
                    
        # discard z coord
        relative_vectors = (self.buffer_coords[tau_index][oxygens_filter][:,:2] - self.buffer_coords[0][oxygens_filter][:,:2])
        
        squared_disp = np.linalg.norm(relative_vectors, axis=1)**2
        
        # bin_index, tau_index, sd, n_entries
    
        return bin_index, tau_index, np.sum(squared_disp), squared_disp.shape[0]
        
    def finalise_analysis(self, save=True, filename="binned_msd"):

        # nan_filter = (~(np.isnan(bin_msds)) & ~(np.isinf(bin_msds)) & ())
        
        bin_msds = self.cum_bin_msd/self.cum_bin_entries
        
        taus = np.arange(self.tau_range, step=self.tau_step)    
        
        ms = np.array([])
        cs = np.array([])
        
        straight_line = lambda x, m: m*x
        
        for msd in bin_msds:
            
            #aaaaaaaaaaaaaaaaa this doesn't work fix it later
            
            
            nan_filter = (~(np.isnan(msd)) & ~(np.isinf(msd)) & ~(np.isclose(msd, 0)))
            
            if not nan_filter.any():
                ms = np.append(ms, 0)
                # cs = np.append(cs, 0)
                continue
            
            timesteps = np.array(taus)[nan_filter]
            msd = msd[nan_filter]
            
            
            popt, pcov = curve_fit(straight_line, timesteps, msd)
            
            print(popt)
            
            ms = np.append(ms, popt[0])
            # cs = np.append(cs, popt[1])
            
        # ms_filter = ms != 0
        plot_bins = self.bins # [ms_filter]
        plot_ms = self.bins # [ms_filter]
        Ds = (plot_ms/4)*1e-5
        np.savetxt(f"{filename}.txt", (plot_bins, Ds))

        plt.plot(plot_bins, Ds)
        plt.xlabel("Z Coordinate (Å)")
        plt.ylabel("Diffusion coefficient ($m^2/s^{-1}$)")
        plt.tight_layout()
        plt.savefig(f"{filename}.png")
        plt.cla()
        plt.clf()

        # for bin_index, msd in zip(range(self.nbins), bin_msds):

        #     plt.plot(taus, msd)
        #     plt.xlabel("Tau (fs)")
        #     plt.ylabel("MSD Å^2")
        #     plt.tight_layout()
        #     plt.savefig(f"msd_plots/{bin_index}_bin_msd_plot.png")
            # plt.cla()
            # plt.clf()
        
        
class TetrahedralAngleAnalysis(Analysis):
    
    O_types = ["OS"]
    C_types = ["CT"]
    C_O_cutoff=17
    # it's easier to just hardcode them in as they never change
    __q_indices = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    
    def __init__(self, bin_height=0.1, fixed_bins=True, fixed_dims=True):
        super().__init__(bin_height, fixed_bins)
        
        self.O_atoms = None
        self.O_ids = None
        self.n_Os = 0

        self.C_atoms = None
        self.C_ids = None
        self.n_Cs = 0

        self.initialised = False
        
        self.bin_height = bin_height

        self.nbins = int(np.ceil(np.abs(self.C_O_cutoff/self.bin_height)))+1
        self.bin_height = np.abs(self.C_O_cutoff)/(self.nbins-1)


        padding = 20
        self.cum_q_binned = np.zeros((self.nbins+padding,))
        self.n_bin_entries  = np.zeros((self.nbins+padding,))

        self.bins = np.linspace(0, self.C_O_cutoff, self.nbins+padding)
        
        self.fixed_dims = fixed_dims
        self.dims = None


    def initialise(self, new_atoms, dump):

        if not self.initialised:

            self.set_dims(dump)

            self.O_atoms = np.array(list(filter(lambda a: a.type_id in self.O_types, dump.atoms)))
            self.n_Os = len(self.O_atoms)
            # self.O_atoms = {o.atom_id : o for o in self.O_atoms}
            self.O_ids = np.array([o.atom_id for o in self.O_atoms])

            self.C_atoms = np.array(list(filter(lambda a: a.type_id in self.C_types, dump.atoms)))
            self.n_Cs = len(self.C_atoms)
            # self.O_atoms = {o.atom_id : o for o in self.O_atoms}
            self.C_ids = np.array([c.atom_id for c in self.C_atoms])
            self.initialised = True
        
        return True

    def set_dims(self, dump):
        self.dims = np.array([
                                    dump.xbounds[1]-dump.xbounds[0],
                                    dump.ybounds[1]-dump.ybounds[0],
                                    dump.zbounds[1]-dump.zbounds[0],
                                ])
        return self.dims
      
        
    def __call__(self, old_atoms, new_atoms, *, dump):
        
        self.initialise(new_atoms, dump)

        if not self.fixed_dims:
            self.set_dims(dump)

        # Step specific stuff
        
        # TODO: write a container object for coordinates (subclassing np.array maybe?)
        # to save needing to do this each time
        self.O_coords = np.array([a.coords for a in self.O_atoms])
        self.C_coords = np.array([a.coords for a in self.C_atoms])

        for i, central_atom in enumerate(self.O_atoms):
            
            # C_O_distances = np.linalg.norm(C_coords - central_atom.coords, axis=1)
            C_O_distances = self.__pairwise_distances_pbc(central_atom.coords, self.C_coords)

            assert len(C_O_distances) == self.n_Cs

            # This is only relevant for whether to count the atom/what bin it goes in,
            # so we don't care about which carbon atom it is
            min_distance = np.min(C_O_distances)
            
            if min_distance > self.C_O_cutoff:
                continue

            q = self.get_q(i)

            bin_index = np.floor(min_distance/self.bin_height).astype(int)

            self.n_bin_entries[bin_index] += 1
            self.cum_q_binned[bin_index] += q

    def finalise_analysis(self, save=True, filename="binned_tetrahedral"):

        # x==y
        mean_q = (self.cum_q_binned[:len(self.bins)]/self.n_bin_entries[:len(self.bins)])
        np.savetxt(f"{filename}.txt", (self.bins, mean_q))

        datestring = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")

        plt.plot(self.bins, mean_q)
        plt.xlabel("Distance from central carbon (Å)")
        plt.ylabel("Tetrahedral order parameter")
        plt.title(f"Generated at {datestring}")
        plt.tight_layout()
        plt.savefig(f"{filename}.png")
        plt.cla()
        plt.clf()

        print("thanks for playing :+)")

    def __pairwise_distances_pbc(self, central_pos, pos_list):

        # For subtracting the dims in an efficient way with numpy it's easier to 
        # have things in the same shape, although this makes it a little less efficient
        tiled_dims = np.tile(self.dims, (len(pos_list), 1))

        # given we're just trying to calculate an absolute distance we don't care about the sign
        # and only the magnitude. Therefore we can safely ignore that, and only work with the
        # absolute value
        relative_vectors = np.abs(pos_list - central_pos)
        relative_vectors[relative_vectors > self.dims/2] -= tiled_dims[relative_vectors > self.dims/2]

        return np.linalg.norm(relative_vectors, axis=1)
    
    def __relative_vectors_pbc(self, central_pos, nearest_neighbours):
        
        pbc_vectors = nearest_neighbours - central_pos
        tiled_dims = np.tile(self.dims, (len(nearest_neighbours), 1))

        pbc_vectors[pbc_vectors > self.dims/2] -= tiled_dims[pbc_vectors > self.dims/2]
        pbc_vectors[pbc_vectors < -self.dims/2] += tiled_dims[pbc_vectors < -self.dims/2]

        return pbc_vectors

    def get_q(self, central_index):

        central_O_coords = self.O_coords[central_index]
        distances = self.__pairwise_distances_pbc(central_O_coords, self.O_coords)
            
        assert len(distances) == self.n_Os

        # Doing it this way as the lowest difference will always be 0 for i == j
        nearest_indices = np.argpartition(distances, 5)[:5]
        # Remove own atom
        nearest_indices = nearest_indices[nearest_indices != central_index]

        # half assed inline unit testing, disable with -O flag
        assert len(nearest_indices) == 4
        
        nearest_neighbours = self.O_atoms[nearest_indices]
        nearest_neighbours_coords = np.array(list(map(lambda x: x.coords, nearest_neighbours)))
        
        #vecs = np.array(nearest_neighbours) - central_O
        vecs = self.__relative_vectors_pbc(central_O_coords, nearest_neighbours_coords)
        
        # pairs = list(zip(vecs, np.roll(vecs, -1)))[:-1]
        
        temp_cum = 0
        
        for i, j in self.__q_indices:

            vec1 = vecs[i]
            vec2 = vecs[j]

            cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))

            temp_cum += (cos_theta + 1/3)**2
        
        return 1 - (3/8)*temp_cum
    

class TetrahedralSolvationShellAnalysis(TetrahedralAngleAnalysis):
    
    C_types = np.array(["CT", "CC", "CT2"])
    shell_boundaries = np.array([5, 10, 15])
    
    def __init__(self, bin_height=0.1, fixed_bins=True, fixed_dims=True):
        super().__init__(bin_height, fixed_bins, fixed_dims)
        
        self.C_O_cutoff = self.shell_boundaries[0]
        
        self.shell_cum = np.zeros(self.shell_boundaries.shape)
        self.shell_n_entries = np.zeros(self.shell_boundaries.shape)
        
    def __call__(self, old_atoms, new_atoms, *, dump):
        
        self.initialise(new_atoms, dump)

        if not self.fixed_dims:
            self.set_dims(dump)
            

        # Step specific stuff
        
        # TODO: write a container object for coordinates (subclassing np.array maybe?)
        # to save needing to do this each time
        self.O_coords = np.array([a.coords for a in self.O_atoms])
        self.C_coords = np.array([a.coords for a in self.C_atoms])

        for i, central_atom in enumerate(self.O_atoms):
            
            C_O_distances = self.__pairwise_distances_pbc(central_atom.coords, self.C_coords)

            assert len(C_O_distances) == self.n_Cs

            # This is only relevant for whether to count the atom/what bin it goes in,
            # so we don't care about which carbon atom it is
            min_distance = np.min(C_O_distances)
            
            if min_distance > self.C_O_cutoff:
                continue

            q = self.get_q(i)

            shell_index = np.where(min_distance <= self.shell_boundaries)[0][0]
            
            self.shell_cum[shell_index] += q 
            self.shell_n_entries[shell_index] += 1

    def finalise_analysis(self, save=True, filename="tetrahedral_distribution"):

        mean_shell_q = self.shell_cum/self.shell_n_entries
        
        np.savetxt(f"{filename}.txt", (self.shell_boundaries, np.arange(0, len(self.shell_boundaries)), mean_shell_q))

        datestring = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")

        print(f"Completed at {datestring}")
        for i, mean_q, shell_hi in enumerate(zip(mean_shell_q, self.shell_boundaries))):
            
            if i == 0:
                shell_lo = 0
            else:
                shell_lo = self.shell_boundaries[i-1]
            
            print(f"Solvation shell {i + 1} ({shell_lo} -> {shell_hi} Å) mean q: {mean_q}")
            
            
        print("thanks for playing :+)")
        

class TetrahedralDistributionAnalysis(TetrahedralAngleAnalysis):

    def __init__(self, bin_height=0.01, fixed_bins=True, fixed_dims=True, outfile=None):
        super().__init__(bin_height, fixed_bins, fixed_dims)

        self.nbins = 1/self.bin_height
        assert int(self.nbins) == self.nbins
        self.nbins = int(self.nbins)

        self.bins = np.linspace(0, 1, self.nbins)
        self.n_bin_entries = np.zeros((self.nbins,))
        self.outfile=outfile

    def __call__(self, old_atoms, new_atoms, *, dump):
        
        self.initialise(new_atoms, dump)

        if not self.fixed_dims:
            self.set_dims(dump)
            
        if self.outfile is not None and dump.current_step % 2000000 == 0:
            self.finalise_analysis(filename=self.outfile)

        # Step specific stuff
        
        # TODO: write a container object for coordinates (subclassing np.array maybe?)
        # to save needing to do this each time
        self.O_coords = np.array([a.coords for a in self.O_atoms])
        self.C_coords = np.array([a.coords for a in self.C_atoms])

        for i, central_atom in enumerate(self.O_atoms):

            q = self.get_q(i)

            bin_index = np.floor(q/self.bin_height).astype(int)

            self.n_bin_entries[bin_index] += 1

    def finalise_analysis(self, save=True, filename="tetrahedral_distribution"):


        np.savetxt(f"{filename}.txt", (self.bins, self.n_bin_entries))

        datestring = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")

        plt.plot(self.bins, self.n_bin_entries)
        plt.xlabel("Tetrahedral Angle Order Parameter")
        plt.ylabel("Relative proportion (arbitrary units)")
        plt.title(f"Generated at {datestring}")
        plt.tight_layout()
        plt.savefig(f"{filename}.png")
        plt.cla()
        plt.clf()

        print("thanks for playing :+)")