import numpy as np
from numpy import ndarray
import pandas as pd
from scipy.optimize import linear_sum_assignment
from pathlib import Path
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import cpu_count
from functools import partial

atomic_masses = {
    'Mo': 95.95,
    'Cl': 35.45,
    'C': 12.01,
    'S': 32.07,
    'F': 18.998,
    'B': 10.81,
    'Ag': 107.87,
    'P': 30.97,
    'N': 14.01,
    'Rh': 102.91,
    'Fe': 55.85,
    'W': 183.84,
    'Br': 79.90,
    'O': 16.00,
    'Pd': 106.42,
    'Na': 22.99,
    'H': 1.008
}

atomic_radii = {'Mo': 1.54,
 'Cl': 1.02,
 'C': 0.7683762219460135,
 'S': 1.0761039450679384,
 'F': 0.57,
 'B': 0.8974061360281438,
 'Au': 1.4367724089623388,
 'Ag': 1.475936900884955,
 'P': 0.9831384041006173,
 'N': 0.7366888871898063,
 'Rh': 1.2365903179673157,
 'Fe': 1.3277038874098017,
 'W': 1.62,
 'Br': 1.2,
 'O': 0.638258821102083,
 'Pd': 1.39,
 'Na': 1.66,
 'H': 0.372}

#Alpha from PDF script
alpha = 1.25

# Constants for conversion and calculations
h = 6.62607015e-34              # Planck constant in Js
c = 2.99792458e10               # Speed of light in cm/s
Da_to_kg = 1.66053906660e-27    # Conversion factor from Da to kg
A_to_m = 1e-10                  # Conversion factor from A to m
pi = np.pi                      # Value of pi

def parse_xyz_to_df(file_path):
    """
    Parse an XYZ file and return the number of atoms and a DataFrame with atomic coordinates.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()[2:]             # Skip the first two lines
        data = []                             
        for line in lines:
            parts = line.split()              
            atom = parts[0]                   
            x, y, z = map(float, parts[1:4])  
            data.append([atom, x, y, z])      
    
    df = pd.DataFrame(data, columns=['Atom', 'X', 'Y', 'Z'])
    num_atoms = len(df)
    return num_atoms, df

def compute_center_of_mass(df):
    """
    Compute the weighted center of mass of the molecule.
    """
    total_mass = 0     
    com = np.zeros(3)  
    for _, row in df.iterrows():
        mass = atomic_masses[row['Atom']]   
        com += mass * row[['X', 'Y', 'Z']]  
        total_mass += mass                  
    com /= total_mass                       
    return com  

def translate_to_com(df):
    """
    Translate the molecule to its center of mass.
    """
    com = compute_center_of_mass(df)
    df_translated = df.copy()
    df_translated[['X', 'Y', 'Z']] -= com
    return df_translated

def inertia_moment_eigenval_eigenvec(df):
    """
    Calculate the eigenvalues and eigenvectors of the inertia moment tensor.
    """
    X = df[['X', 'Y', 'Z']].values.astype(np.float64)
    I_matrix = np.zeros((3, 3), dtype=np.float64)
    
    for i, atom in enumerate(df['Atom']):
        mass = atomic_masses[atom]         
        I_matrix += mass * np.outer(X[i], X[i])
    
    eigenvalues, eigenvectors = np.linalg.eigh(I_matrix)
    return eigenvalues, eigenvectors

def rotate_with_eigenvectors(df, eigenvectors):
    """
    Rotate the molecule using the provided eigenvectors.
    """
    rotated_coords = np.dot(df[['X', 'Y', 'Z']].values, eigenvectors)
    df_rotated = df.copy()
    df_rotated[['X', 'Y', 'Z']] = rotated_coords
    return df_rotated

def mirror_image(coords):
    return -coords


def create_df_from_atoms(atom_types, coordinates):
    """
    Create a DataFrame from the given number of atoms, atom types, and coordinates.
    """
    data = {'Atom': atom_types, 'X': coordinates[:, 0], 'Y': coordinates[:, 1], 'Z': coordinates[:, 2]}
    df = pd.DataFrame(data)
    return df

def calculate_neighbors(df, alpha):
    """
    Calculate neighbors for atomic environment analysis.
    """
    df_copy = df.copy()
    df_copy['Radius'] = df_copy['Atom'].map(atomic_radii)
    
    distances = pd.DataFrame(
        np.sqrt(
            (df_copy['X'].values[:, np.newaxis] - df_copy['X'].values)**2 + 
            (df_copy['Y'].values[:, np.newaxis] - df_copy['Y'].values)**2 +
            (df_copy['Z'].values[:, np.newaxis] - df_copy['Z'].values)**2
        )
    )
    
    cutoff_distances = alpha * (
        df_copy['Radius'].values[:, np.newaxis] + df_copy['Radius'].values
    )
    
    neighbors = (distances < cutoff_distances).astype(int)
    np.fill_diagonal(neighbors.values, 0)
    neighbor_list = [list(np.where(row)[0]) for row in neighbors.values]
    
    return neighbor_list

def create_sorted_atom_neighbor_strings(atom_types, neighboring_list):
    """
    Creates a list of strings where each string starts with an atom type followed by its sorted neighboring atom types.
    """
    atom_neighbor_strings = []
    for atom_index, atom_type in enumerate(atom_types):
        neighbor_types = [atom_types[i] for i in neighboring_list[atom_index]]
        sorted_neighbor_types = sorted(neighbor_types)
        atom_neighbor_strings.append(atom_type + ''.join(sorted_neighbor_types))
    return atom_neighbor_strings

def hausdorff_minimal_distance_by_type(P, Q, atom_types_P, atom_types_Q):
    """
    Calculate Hausdorff distance while considering minimal distance and atom types using optimal atom assignment.
    """
    unique_atom_types = set(atom_types_P)
    max_hausdorff_distance = 0  

    for atom_type in unique_atom_types:
        indices_P = [i for i, at in enumerate(atom_types_P) if at == atom_type]
        indices_Q = [i for i, at in enumerate(atom_types_Q) if at == atom_type]

        if len(indices_P) == 0 or len(indices_Q) == 0:
            continue  

        P_subset = P[indices_P]
        Q_subset = Q[indices_Q]

        from scipy.spatial.distance import cdist
        distances = cdist(P_subset, Q_subset)
        min_dist_P_to_Q = np.min(distances, axis=1)  
        min_dist_Q_to_P = np.min(distances, axis=0)  
        atom_type_hausdorff_distance = max(min_dist_P_to_Q.max(), min_dist_Q_to_P.max())
        max_hausdorff_distance = max(max_hausdorff_distance, atom_type_hausdorff_distance)

    return max_hausdorff_distance


def rmsd_minimal_distance_by_type(P, Q, atom_types_P, atom_types_Q):
    """
    Calculate RMSD (Root Mean Square Deviation) while considering minimal distance and atom types using optimal atom assignment.
    
    :param P: numpy.ndarray, a set of points representing one molecular structure (e.g., atomic coordinates).
    :param Q: numpy.ndarray, a set of points representing another molecular structure to be compared with P.
    :param atom_types_P: list or numpy.ndarray, atom types corresponding to the points in P.
    :param atom_types_Q: list or numpy.ndarray, atom types corresponding to the points in Q.
    :return: float, the calculated RMSD value based on the optimal atom assignment considering atom types.
    """
    unique_atom_types = set(atom_types_P)  # Assuming atom_types_P and atom_types_Q have the same unique elements
    total_squared_distance = 0
    total_atoms = 0

    for atom_type in unique_atom_types:
        indices_P = [i for i, at in enumerate(atom_types_P) if at == atom_type]
        indices_Q = [i for i, at in enumerate(atom_types_Q) if at == atom_type]

        if len(indices_P) == 0 or len(indices_Q) == 0:
            continue  # Skip if no atoms of this type are present in either molecule

        P_subset = P[indices_P]
        Q_subset = Q[indices_Q]

        # Calculate pairwise Euclidean distances between each pair of points of the same atom type in P and Q.
        distances = np.sqrt(np.sum((P_subset[:, np.newaxis, :] - Q_subset[np.newaxis, :, :]) ** 2, axis=2))

        # Use the Hungarian algorithm to find the optimal assignment between points of the same atom type in P and Q.
        row_ind, col_ind = linear_sum_assignment(distances)

        # Accumulate the squared distances based on the optimal assignment of points.
        total_squared_distance += np.sum((P_subset[row_ind] - Q_subset[col_ind]) ** 2)
        total_atoms += len(row_ind)

    if total_atoms > 0:
        rmsd = np.sqrt(total_squared_distance / total_atoms)
    else:
        rmsd = np.nan  # Return NaN if there are no atoms to compare

    # Return the calculated RMSD.
    return rmsd

def molecule_diameter(coordinates):
    """
    Calculate the diameter of a molecule (maximum interatomic distance).
    """
    n_atoms = coordinates.shape[0]
    max_distances = np.zeros(n_atoms)

    for i in range(n_atoms):
        diffs = coordinates - coordinates[i]
        distances = np.linalg.norm(diffs, axis=1)
        max_distances[i] = np.max(distances)
    
    diameter = np.max(max_distances)
    print('Molecule diameter', diameter)
    return diameter

def compute_rotation_batch_with_env(args):
    """
    Compute RMSD for a batch of rotation angles using atomic environment information.
    Uses standard Z-Y-X (yaw-pitch-roll) rotation convention.
    """
    original_coords, mirror_coords, atom_neighbor_strings, atom_neighbor_strings_mirror, angle_batch = args
    
    batch_size = len(angle_batch)
    min_rmsd = np.inf
    optimal_rotation = None
    optimal_coords = None
    
    # Convert angles to arrays for vectorized computation
    angles_array = np.array(angle_batch)
    angle_x_batch = angles_array[:, 0]  # roll
    angle_y_batch = angles_array[:, 1]  # pitch
    angle_z_batch = angles_array[:, 2]  # yaw
    
    # Vectorized trigonometric calculations
    cx, cy, cz = np.cos(angle_x_batch), np.cos(angle_y_batch), np.cos(angle_z_batch)
    sx, sy, sz = np.sin(angle_x_batch), np.sin(angle_y_batch), np.sin(angle_z_batch)
    
    # Process in smaller sub-batches to manage memory
    sub_batch_size = min(100, batch_size)
    
    for i in range(0, batch_size, sub_batch_size):
        end_idx = min(i + sub_batch_size, batch_size)
        
        # Get sub-batch indices
        cx_sub = cx[i:end_idx]
        cy_sub = cy[i:end_idx]
        cz_sub = cz[i:end_idx]
        sx_sub = sx[i:end_idx]
        sy_sub = sy[i:end_idx]
        sz_sub = sz[i:end_idx]
        
        # Vectorized rotation matrix computation (Z-Y-X)
        n_sub = end_idx - i
        rotation_matrices = np.zeros((n_sub, 3, 3))
        
        rotation_matrices[:, 0, 0] = cz_sub * cy_sub
        rotation_matrices[:, 0, 1] = cz_sub * sy_sub * sx_sub - sz_sub * cx_sub
        rotation_matrices[:, 0, 2] = cz_sub * sy_sub * cx_sub + sz_sub * sx_sub
        
        rotation_matrices[:, 1, 0] = sz_sub * cy_sub
        rotation_matrices[:, 1, 1] = sz_sub * sy_sub * sx_sub + cz_sub * cx_sub
        rotation_matrices[:, 1, 2] = sz_sub * sy_sub * cx_sub - cz_sub * sx_sub
        
        rotation_matrices[:, 2, 0] = -sy_sub
        rotation_matrices[:, 2, 1] = cy_sub * sx_sub
        rotation_matrices[:, 2, 2] = cy_sub * cx_sub
        
        # Apply rotations to all coordinates at once
        rotated_coords = np.einsum('ijk,mk->imj', rotation_matrices, mirror_coords)
        
        # Compute RMSD for each rotation in the sub-batch using atomic environment
        for j in range(n_sub):
            current_rmsd = rmsd_minimal_distance_by_type(
                original_coords, 
                rotated_coords[j], 
                atom_neighbor_strings, 
                atom_neighbor_strings_mirror
            )
            
            if current_rmsd < min_rmsd:
                min_rmsd = current_rmsd
                optimal_rotation = angle_batch[i + j]
                optimal_coords = rotated_coords[j]
    
    return min_rmsd, optimal_rotation, optimal_coords


def compute_optimal_rotation_and_rmsd_parallel(original_coords, mirror_coords, atom_neighbor_strings, atom_neighbor_strings_mirror, angle_increment, n_processes=None):
    """
    Parallelized version of compute_optimal_rotation_and_rmsd with atomic environment sensitivity.
    """
    if n_processes is None:
        n_processes = cpu_count()
    
    angle_increment_rad = np.radians(angle_increment)
    angles_360 = np.arange(0, 2*np.pi, angle_increment_rad)
    angles_180 = np.arange(0, np.pi, angle_increment_rad)
    
    print(f"Exploring {len(angles_360) * len(angles_180) * len(angles_360)} angle combinations")
    print(f"Using {n_processes} processes with atomic environment sensitivity")
    
    # Generate all angle combinations
    all_angles = []
    for theta_x in angles_360:  # roll
        for theta_y in angles_180:  # pitch
            for theta_z in angles_360:  # yaw
                all_angles.append((theta_x, theta_y, theta_z))
    
    # Split into larger batches for process-level parallelization
    batch_size = max(1000, len(all_angles) // n_processes)
    angle_batches = [all_angles[i:i + batch_size] for i in range(0, len(all_angles), batch_size)]
    
    print(f"Created {len(angle_batches)} batches of ~{batch_size} angles each")
    
    # Prepare arguments for each process - now includes atomic environment data
    args_list = [(original_coords, mirror_coords, atom_neighbor_strings, atom_neighbor_strings_mirror, batch) 
                 for batch in angle_batches]
    
    # Process batches in parallel using processes
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        results = list(executor.map(compute_rotation_batch_with_env, args_list))
    
    # Find the best result across all batches
    best_result = min(results, key=lambda x: x[0])
    print(f"Best RMSD found: {best_result[0]:.6f}")
    
    return best_result

def process_single_file_with_atomic_env(geom_path, angle_step, alpha_value, n_processes=None):
    """
    Process a single XYZ file with atomic environment sensitivity and parallelization.
    """
    if n_processes is None:
        n_processes = cpu_count()
    
    print(f"Processing {geom_path} with {n_processes} processes and alpha = {alpha_value}")
    
    try:
        natms, df_MOL = parse_xyz_to_df(geom_path)
        df_MOL_CM = translate_to_com(df_MOL)
        atom_types = df_MOL['Atom'].tolist()
        
        # Calculation of the inertia moment and alignment with respect to the principal axes
        eigenval, eigenvec = inertia_moment_eigenval_eigenvec(df_MOL_CM)
        df_MOL_aling_eienvec = rotate_with_eigenvectors(df_MOL_CM, eigenvec)
        MOL_aling_eienvec = np.array(df_MOL_aling_eienvec[['X', 'Y', 'Z']], dtype=np.float64)
        
        # Create mirror image
        MOL_mirror = mirror_image(MOL_aling_eienvec)
        df_MOL_mirror = create_df_from_atoms(atom_types, MOL_mirror)
        
        # Calculate atomic environments for both original and mirror structures
        neighbors = calculate_neighbors(df_MOL, alpha_value)
        neighbors_mirror = calculate_neighbors(df_MOL_mirror, alpha_value)
        atom_neighbor_strings = create_sorted_atom_neighbor_strings(atom_types, neighbors)
        atom_neighbor_strings_mirror = create_sorted_atom_neighbor_strings(atom_types, neighbors_mirror)
        
        # Parallelized optimization with atomic environment sensitivity
        min_rmsd, optimal_rotation, optimal_rotated_coords = compute_optimal_rotation_and_rmsd_parallel(
            original_coords=MOL_aling_eienvec,
            mirror_coords=MOL_mirror,
            atom_neighbor_strings=atom_neighbor_strings,
            atom_neighbor_strings_mirror=atom_neighbor_strings_mirror,
            angle_increment=angle_step,
            n_processes=n_processes
        )

        # If needed the descriptor can be normalized by molecule diameter
        # normalized = min_rmsd / molecule_diameter(MOL_aling_eienvec)
        return min_rmsd
        
    except Exception as e:
        print(f'Error processing {geom_path}: {e}')
        return None

def HCM_atomic_env_parallel(path, angle_step, alpha_value=1.25, n_processes=None):
    """
    Parallelized HCM calculation with atomic environment sensitivity.
    Combines the best of both approaches: atomic environment analysis + parallelization.
    """
    if n_processes is None:
        n_processes = cpu_count()
    
    os.chdir(path)
    xyz_files = list(Path('.').glob('*.xyz'))
    print(f"Found {len(xyz_files)} XYZ files")
    print(f"Using atomic environment sensitivity with alpha = {alpha_value}")
    
    dict_mrmsd = {}
    
    if len(xyz_files) == 1:
        # Single file: use all processes for angle exploration
        geom_path = str(xyz_files[0].resolve())
        key = xyz_files[0].name
        normalized = process_single_file_with_atomic_env(geom_path, angle_step, alpha_value, n_processes)
        if normalized is not None:
            dict_mrmsd[key] = normalized
            print(f'Done {key}: {normalized:.6f}')
    else:
        # Multiple files: balance between file parallelization and angle exploration
        processes_per_file = max(1, n_processes // min(len(xyz_files), 4))
        file_processes = min(n_processes, len(xyz_files))
        
        print(f"Using {file_processes} file processes, {processes_per_file} processes per file")
        
        geom_paths = [str(path.resolve()) for path in xyz_files]
        
        # Create partial function with fixed parameters
        process_func = partial(process_single_file_with_atomic_env, 
                             angle_step=angle_step, 
                             alpha_value=alpha_value,
                             n_processes=processes_per_file)
        
        # Process files in parallel
        with ProcessPoolExecutor(max_workers=file_processes) as executor:
            results = executor.map(process_func, geom_paths)
            
            for i, normalized in enumerate(results):
                if normalized is not None:
                    key = xyz_files[i].name
                    dict_mrmsd[key] = normalized
                    print(f'Done {key}: {normalized:.6f}')
    
    return dict_mrmsd
