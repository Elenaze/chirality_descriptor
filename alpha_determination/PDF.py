import numpy as np
import pandas as pd
from datetime import datetime   
from scipy import stats, signal
from pathlib import Path
from typing import Dict, List, Tuple, Set
from results import AnalysisResults
from matplotlib import pyplot as plt
import seaborn as sns

class PDFAnalyzer:
    def __init__(self, xyz_file: str):
        self.xyz_file = Path(xyz_file)
        self.systems: Dict[int, pd.DataFrame] = {}
        self.unique_atoms: Set[str] = set()
        self.pair_distances: Dict[Tuple[str, str], List[float]] = {}
        self.kde_data: Dict[Tuple[str, str], pd.DataFrame] = {}
        self.discarded_pairs: List[Tuple[str, str]] = []
        self.X_max: Dict[Tuple[str, str], float] = {} 
        self.X_min: Dict[Tuple[str, str], float] = {}
        
        # Reference radii from ASE
        self.ASE_radii = {'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be':0.96, 'B':0.84, 'C':0.76, 'N': 0.71, 
                          'O': 0.66, 'F': 0.57, 'Ne':0.58, 'Na':1.66,'Mg': 1.41, 'Al':1.21, 'Si':1.11, 
                          'P': 1.07,  'S':1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 
                          'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 
                          'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 
                          'Kr': 1.16, 'Rb': 2.20,'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo':1.54, 'Tc':1.47, 
                          'Ru': 1.46,'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In':1.42, 'Sn': 1.39, 'Sb':1.39, 
                          'Te':1.38, 'I': 1.39, 'Xe' :1.40, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce':2.04,'Pr':2.03, 'Nd': 2.01, 
                          'Pm': 1.99, 'Sm':1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb':1.94,'Dy': 1.92,'Ho':1.92, 'Er': 1.89, 'Tm':1.90, 
                          'Yb':1.87,'Lu':1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,'Re': 1.51, 'Os':1.44, 'Ir': 1.41,'Pt': 1.36, 
                          'Au':1.36,'Hg':1.32, 'Tl':1.45, 'Pb':1.46,'Bi':1.48,'Po': 1.40, 'At': 1.50, 'Rn':1.50, 'Fr': 2.60, 'Ra':2.21, 
                          'Ac':2.15, 'Th':2.06, 'Pa': 2.00, 'U':1.96, 'Np': 1.90, 'Pu':1.87,'Am': 1.80, 'Cm':1.69}
        
    def load_systems(self) -> None:
        """Load all systems from XYZ file into DataFrames"""
        natoms = []
        
        # First get number of atoms for each system
        with open(self.xyz_file, 'r') as f:
            for line in f:
                if line.strip().isdigit():
                    natoms.append(int(line.strip()))
        
        # Now extract coordinates
        start = 2
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()
            for n, atoms in enumerate(natoms[::2]):
                end = start + atoms
                data = []
                
                for line in lines[start:end]:
                    atom, x, y, z = line.split()
                    data.append([atom, float(x), float(y), float(z)])
                
                self.systems[n] = pd.DataFrame(
                    data, 
                    columns=['Atom type', 'x', 'y', 'z']
                )
                start = end + 2
    
    def find_unique_atoms(self) -> Set[str]:
        """Get set of unique atoms across all systems"""
        for system in self.systems.values():
            self.unique_atoms.update(system['Atom type'].unique())
        return self.unique_atoms
    
    def calculate_all_pair_distances(self) -> Dict[Tuple[str, str], List[float]]:
        """Calculate distances for all possible atom type pairs"""
        atoms_list = list(self.unique_atoms)
        
        for i, atom1 in enumerate(atoms_list):
            for j in range(i, len(atoms_list)):  # Include self-pairs
                atom2 = atoms_list[j]
                distances = []
                
                # Loop through all systems
                for system in self.systems.values():
                    coords1 = system[system['Atom type'] == atom1][['x','y','z']].values
                    coords2 = system[system['Atom type'] == atom2][['x','y','z']].values
                    
                    # Calculate distances between all pairs
                    for c1 in coords1:
                        for c2 in coords2:
                            if not np.array_equal(c1, c2):  # Skip same atom
                                dist = np.linalg.norm(c1 - c2)
                                distances.append(dist)
                
                if distances:  # Store only if pairs were found
                    self.pair_distances[(atom1, atom2)] = distances
        
        return self.pair_distances
    
    def compute_kde(self, bw_method: float = 0.04) -> None:
        """Compute KDE for all atom pairs"""
        for key, distances in self.pair_distances.items():
            kde = stats.gaussian_kde(distances, bw_method=bw_method)
            x = np.linspace(min(distances), max(distances), 1000)
            y = kde(x)
            self.kde_data[key] = pd.DataFrame({'X': x, 'Y': y})

    
    def find_max_min(self, min_peak_height: float = 0.002) -> Tuple[Dict[Tuple[str, str], float], Dict[Tuple[str, str], float]]:
        """Find first significant maximum and following minimum in KDE for each atom pair"""
        for key in self.kde_data:
            # Find maxima
            max_index = signal.find_peaks_cwt(self.kde_data[key]['Y'], np.arange(1,10))
            if max_index.size != 0:
                for idx in max_index:
                    if self.kde_data[key].iloc[idx]['Y'] > min_peak_height:
                        x_max = self.kde_data[key].iloc[idx]['X']
                        self.X_max[key] = x_max
                        break
            
            # Find minima after maximum
            if key in self.X_max:
                min_index, _ = signal.find_peaks(-self.kde_data[key]['Y'])
                if min_index.size != 0:
                    for idx in min_index:
                        if self.kde_data[key].iloc[idx]['X'] > self.X_max[key]:
                            x_min = self.kde_data[key].iloc[idx]['X']
                            self.X_min[key] = x_min
                            break

        return self.X_max, self.X_min
    

    def determine_alpha(self, radii: Dict[str, float],
                       min_samples: int = 200,
                       max_alpha: float = 2.0) -> pd.DataFrame:
        """Calculate alpha values using stored maxima and minima"""
        alpha_values = {}

        for key in self.kde_data:
            if len(self.pair_distances[key]) < min_samples:
                print(f"Not enough statistics for {key}")
                continue

            if key not in self.X_min or key not in self.X_max:
                continue

            alpha = self.X_min[key] / (radii[key[0]] + radii[key[1]])
            
            if alpha <= max_alpha:
                if alpha <= (radii[key[0]] + radii[key[1]]):
                    alpha_values[key] = alpha
                else:
                    alpha_values[key] = np.nan
                    
            print(f"Pair: {key}, Alpha: {alpha:.3f}")

        # Create DataFrame and sort
        df = pd.DataFrame({'Pair': list(alpha_values.keys()), 
                          'Alpha': list(alpha_values.values())})
        df.sort_values(by='Alpha', inplace=True)
        
        return df.dropna()
    

    def determine_same_atom_radii(self, min_atom_count: int = 100, min_peak_height: float = 0.002, max_alpha: float = 1.2) -> Dict[str, List[float]]:
        """Determine atomic radii from same-atom pair distances"""

        radii = {}
        self.discarded_pairs = []
    
        for key in self.kde_data:
            if key[0] != key[1]:  # Skip different atom pairs
                continue
            
            atom_type = key[0]
            atom_count = sum(1 for sys in self.systems.values() 
                        for atom in sys['Atom type'] if atom == atom_type)
        
        # Check statistics
            if atom_count < min_atom_count:
                print(f"Not enough {atom_type} atoms ({atom_count})")
                radii[atom_type] = [0]
                self.discarded_pairs.append(key)
                continue
            
        # Find first maximum
            max_index, _ = signal.find_peaks(self.kde_data[key]['Y'])
            x_max = None
        
            if max_index.size != 0:
                for idx in max_index:
                    if self.kde_data[key].iloc[idx]['Y'] > min_peak_height:
                        x_max = self.kde_data[key].iloc[idx]['X']
                        break
                    
                if x_max is not None:
                    alpha = x_max / (2 * self.ASE_radii[atom_type])
                    if alpha <= max_alpha:
                        radii[atom_type] = [x_max / 2]
                        print(f"Radius for {atom_type}: {x_max/2:.3f} Å "
                              f"(ASE: {self.ASE_radii[atom_type]:.3f} Å)")
                    else:
                        print(f"Distance too large for {atom_type}")
                        radii[atom_type] = [0]
                        self.discarded_pairs.append(key)
                    
        return radii

    def determine_different_atom_radii(self, radii: Dict[str, List[float]],
                                 min_samples: int = 100,
                                 min_peak_height: float = 0.002,
                                 max_alpha: float = 1.1) -> Dict[str, List[float]]:
        """
    Determine remaining atomic radii using different-atom pairs
    
    Args:
        radii: Dictionary of already determined radii
        min_samples: Minimum number of distance samples
        min_peak_height: Minimum peak height for detection
        max_alpha: Maximum allowed ratio vs ASE radii
        
    Returns:
        Updated radii dictionary
        """
        for key in self.kde_data:
            if key[0] == key[1]:  # Skip same atom pairs
                continue
            
            atom1, atom2 = key
            if ((radii[atom1][0] != 0 and radii[atom2][0] == 0) or
                (radii[atom1][0] == 0 and radii[atom2][0] != 0)):
            
            # Check statistics
                if len(self.pair_distances[key]) < min_samples:
                    print(f"Not enough samples for {key}")
                    self.discarded_pairs.append(key)
                    continue
                
            # Find first maximum
                max_index, _ = signal.find_peaks(self.kde_data[key]['Y'])
                if max_index.size == 0:
                    continue
                
                for idx in max_index:
                    if self.kde_data[key].iloc[idx]['Y'] > min_peak_height:
                        x_max = self.kde_data[key].iloc[idx]['X']
                        break
                else:
                    continue
                
            # Calculate new radius
                if radii[atom1][0] == 0:
                    alpha = x_max / (self.ASE_radii[atom1] + radii[atom2][0])
                    if alpha <= max_alpha:
                        new_radius = x_max - radii[atom2][0]
                        radii[atom1].append(new_radius)
                        print(f"Radius for {atom1}: {new_radius:.3f} Å")
                    else:
                        self.discarded_pairs.append(key)
                else:
                    alpha = x_max / (radii[atom1][0] + self.ASE_radii[atom2])
                    if alpha <= max_alpha:
                        new_radius = x_max - radii[atom1][0]
                        radii[atom2].append(new_radius)
                        print(f"Radius for {atom2}: {new_radius:.3f} Å")
                    else:
                        self.discarded_pairs.append(key)
                    
        return radii
    
    
    def save_results(self, output_dir: str):
        """Save analysis results"""
        results = AnalysisResults(
            timestamp=datetime.now().isoformat(),
            xyz_file=str(self.xyz_file),
            pair_distances=self.pair_distances,
            kde_data=self.kde_data,
            X_max=self.X_max,
            radii=self.radii,
            alpha_values=self.alpha_values
        )
        results.save(output_dir)

    def plot_pair_distributions(self, output_dir: str | None = None):
        """
        Plot pair distribution functions with maxima and minima markers
    
        Args:
            output_dir: Optional directory to save plots
        """
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
        plt.rcParams['font.size'] = 16
    
        bins = np.linspace(0, 10, 101)
    
        for key in self.kde_data.keys():
            if key in self.pair_distances:
            # Create figure
                fig, ax = plt.subplots()
            
            # Plot histogram and KDE
                sns.histplot(self.pair_distances[key], bins=bins, 
                        stat='density', alpha=0.4, element='step')
                ax.plot(self.kde_data[key]['X'], self.kde_data[key]['Y'])
            
            # Set plot style
                ax.set_title(f'{key[1]}-{key[0]}', fontsize=16)
                ax.set_xlabel(r'd / $\mathrm{\AA}$')
                ax.set_facecolor('None')
                ax.set_xlim(0, 10)
            
            # Add vertical lines for maxima and minima
                if key in self.X_max:
                    ax.axvline(x=self.X_max[key], linestyle='--', 
                          color='black', label='X max')
                
                if key in self.X_min:
                    ax.axvline(x=self.X_min[key], linestyle='--',
                          color='salmon', label='X min')
            
                ax.legend()
            
            # Save plot if output directory specified
                if output_dir:
                    out_path = Path(output_dir)
                    out_path.mkdir(parents=True, exist_ok=True)
                    plt.savefig(out_path / f'pdf_{key[0]}_{key[1]}.pdf',
                           bbox_inches='tight', dpi=300)
            
                plt.close()

    
    def analyze(self) -> Dict[Tuple[str, str], List[float]]:
        """Run complete analysis workflow"""
        print("Loading systems...")
        self.load_systems()
        
        print("Finding unique atoms...")
        self.find_unique_atoms()
        
        print("Calculating pair distances...")
        self.calculate_all_pair_distances()
        
        print("Computing KDE...")
        self.compute_kde()
        
        print("Finding maximum and minum...")
        self.find_max_min()
        
        return self.pair_distances