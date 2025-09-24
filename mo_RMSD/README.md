# Molecular Chirality Descriptor (mo_RMSD)

This package provides tools for computing chirality descriptors of molecular structures using a modified Root Mean Square Deviation (RMSD) approach with atomic environment sensitivity.

## Overview

The mo_RMSD (molecular RMSD) chirality descriptor quantifies the degree of chirality in molecular structures by:
- Computing the optimal superposition between a molecule and its mirror image
- Using atomic environment analysis to enhance sensitivity to local chemical environments
- Normalizing by molecular diameter to provide scale-independent descriptors
- Employing parallel processing for efficient exploration of rotational space

## Features

- **Atomic Environment Sensitivity**: Uses atomic neighbor analysis with customizable cutoff radii
- **Parallel Processing**: Multi-core optimization for rotational space exploration
- **Incremental Results Saving**: Results saved immediately after each system completion
- **Resume Capability**: Automatically detects and skips completed analyses
- **Flexible Input**: Supports both batch processing and single system analysis
- **Multiple Output Formats**: JSON, CSV, and legacy DAT formats

## Installation

Ensure you have the required Python packages installed:

```bash
pip install numpy pandas scipy pathlib concurrent.futures
```

## Core Algorithm

### 1. **Molecular Preprocessing**
- Parse XYZ coordinate files
- Translate molecules to center of mass
- Align with principal axes of inertia

### 2. **Atomic Environment Analysis**
- Calculate atomic neighbors using distance cutoffs: `d_cutoff = α × (r_i + r_j)`
- Create atomic environment fingerprints
- Default α = 1.25 (optimized from PDF analysis)

### 3. **Mirror Image Generation**
- Generate enantiomer by coordinate inversion
- Apply systematic rotational search (Euler angles)
- Use Hungarian algorithm for optimal atom assignment

### 4. **Chirality Descriptor Calculation**
- Compute environment-sensitive RMSD for each rotation
- Find minimum RMSD across all rotations
- Normalize by molecular diameter: `χ = RMSD_min / d_molecule`

## Supported Atom Types

The following atoms are supported with predefined atomic masses and radii:

**Common Elements:**
- H (Hydrogen): mass=1.008, radius=0.372 Å
- C (Carbon): mass=12.01, radius=0.768 Å  
- N (Nitrogen): mass=14.01, radius=0.737 Å
- O (Oxygen): mass=16.00, radius=0.638 Å

**Extended Set:**
- Transition metals: Mo, Rh, Fe, W, Pd, Ag, Au
- Halogens: F, Cl, Br
- Other elements: B, P, S, Na

## Usage

### Command Line Interface

#### **Batch Processing (Multiple Systems):**
```bash
python main.py <angle_step> [data_directory] [n_processes]
```

Examples:
```bash
python main.py 30                    # Default: ../data, 8 processes
python main.py 30 data 4            # Custom data dir and processes
```

#### **Single System Analysis:**
```bash
python main.py <angle_step> <data_dir> <n_processes> <system_directory>
```

Examples:
```bash
python main.py 30 data 4 data/SR/cAlaAla      # Analyze specific system
python main.py 90 data 2 data/SS/cPhePhe      # Quick analysis
```

### Configuration File Usage

#### **Create Config File:**
```bash
# config.txt
30                       # angle_step (degrees)
data                     # data_directory  
4                        # n_processes
data/SR/cAlaAla         # system_directory (optional, for single system)
```

#### **Run with Config:**
```bash
python main.py < config.txt
```

#### **Quick Tests:**
```bash
# Fast test with large angle steps
echo -e "90\ndata\n2\ndata/SR/cAlaAla" | python main.py

# Detailed analysis with small angle steps  
echo -e "5\ndata\n8\ndata/SR/cPhePhe" | python main.py
```

## Expected Directory Structure

```
data/
├── SR/                     # S,R chirality
│   ├── cAlaAla/           # Cyclo(Ala-Ala)
│   │   └── *.xyz          # Structure files
│   ├── cGlyGly/           # Cyclo(Gly-Gly)
│   │   └── *.xyz
│   └── ...
└── SS/                     # S,S chirality  
    ├── cAlaAla/
    │   └── *.xyz
    └── ...
```

## Output Files

### **Individual System Results:**
- `{chirality}_{system}_step_{angle}.json` - Complete system data
- Progress tracking and error logging

### **Aggregated Results:**
- `chirality_descriptors_step_{angle}.csv` - Master CSV file
- `{angle}_step_env_hd_{chirality}.dat` - Legacy format
- `chirality_descriptors_complete_step_{angle}.json` - Final summary

### **Monitoring Files:**
- `progress_step_{angle}.txt` - Real-time progress log
- `errors_{angle}.log` - Error log (if any failures occur)

## Parameters

### **Angle Step**
- **5°**: High precision, comprehensive search (~186,624 combinations)
- **30°**: Standard analysis (~1,728 combinations)  
- **90°**: Quick testing (~64 combinations)

### **Process Count**
- **Single system**: 2-8 processes recommended
- **Batch processing**: Use available CPU cores
- **Memory consideration**: More processes = higher memory usage

### **Atomic Environment Parameter (α)**
- **Default**: 1.25 (optimized from PDF analysis)
- **Range**: Typically 1.0-2.0
- **Effect**: Higher values include more distant neighbors

## Performance Guidelines

| Angle Step | Time per System | Precision | Use Case |
|------------|----------------|-----------|----------|
| 90°        | ~5 seconds     | Low       | Quick testing |
| 30°        | ~30 seconds    | Medium    | Standard analysis |
| 10°        | ~5 minutes     | High      | Detailed studies |
| 5°         | ~20 minutes    | Very High | Publication quality |

## Example Results

### **Typical Chirality Descriptor Values:**
- **Highly chiral systems**: χ > 0.15
- **Moderately chiral**: 0.05 < χ < 0.15  
- **Weakly chiral**: χ < 0.05
- **Achiral systems**: χ ≈ 0

### **Sample Output:**
```json
{
    "chirality": "SR",
    "system": "cAlaAla", 
    "results": {
        "2.log-geo-final.xyz": 0.167479
    },
    "processing_time_seconds": 4.51
}
```

## Advanced Features

### **Resume Capability**
- Automatically detects completed systems
- Skips reprocessing in non-interactive mode
- Loads existing results for final summary

### **Error Handling**
- Individual system failures don't stop batch processing
- Detailed error logging with timestamps
- Graceful handling of missing files or invalid structures

### **Memory Management**
- Batch processing of angle combinations
- Efficient vectorized operations
- Configurable sub-batch sizes

## Troubleshooting

### **Common Issues:**

1. **Missing atom types**: Add entries to `atomic_masses` and `atomic_radii` dictionaries
2. **Path resolution errors**: Use absolute paths or ensure correct relative path structure
3. **Memory errors**: Reduce number of processes or angle precision
4. **No results found**: Check XYZ file format and atomic coordinates

### **Performance Optimization:**
- Use larger angle steps for initial screening
- Employ fewer processes for single systems
- Consider system complexity when setting precision

## Scientific Background

The chirality descriptor is based on:
- **Cremer-Pople approach**: Systematic conformational analysis
- **Hungarian algorithm**: Optimal atom-to-atom assignment
- **Atomic environment theory**: Local chemical environment sensitivity
- **Scale normalization**: Diameter-normalized metrics for size independence

## References

- Hausdorff distance applications in molecular comparison
- Root Mean Square Deviation (RMSD) in structural biology
- Atomic environment analysis in chemical informatics
- Parallel optimization algorithms for molecular conformations

## License

This software is provided as-is for research purposes.