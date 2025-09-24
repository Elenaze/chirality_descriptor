# Pair Distribution Function Analysis

Tool for analyzing atomic pair distribution functions from XYZ structure files, with a focus on determining alpha values and atomic radii.

## Overview

This tool was developed to analyze pair distribution functions from molecular structure data, specifically focused on dipeptide conformations. It enables systematic analysis of atomic distances and determination of characteristic parameters.

## Features

- **XYZ File Processing**: Load and analyze molecular structure data
- **PDF Analysis**: Calculate pair distribution functions for all atom pairs
- **Parameter Determination**:
  - Atomic radii from PDF maxima
  - Alpha values from PDF minima 
- **Visualization**:
  - PDF plots with marked maxima/minima
  - Atom abundance histograms
  - Comparative protonation state analysis

## Installation

```bash
# Clone the repository
git clone [your-repo-url]
cd Repository

# Install dependencies
pip install numpy pandas scipy matplotlib seaborn
```

## Usage

Basic example:
```python
from alpha_determination.PDF import PDFAnalyzer

# Initialize analyzer
analyzer = PDFAnalyzer("./data/combined_minima_all.xyz")

# Run complete analysis
analyzer.analyze()

# Get results
df_alpha = analyzer.determine_alpha(radii)
analyzer.plot_pair_distributions("results/plots")
```

## Results Structure

```
results/
└── [timestamp]/
    ├── analysis_results.json  # Main results
    ├── kde_*.csv             # KDE data
    └── full_results.pkl      # Complete data
```

## Method Details

### Atomic Radii Determination
1. Analysis of same-atom pair distributions
2. Refinement using different-atom pairs
3. Validation against ASE reference values

### Alpha Value Calculation
1. PDF maxima/minima identification
2. Ratio calculation: minimum position / sum of atomic radii
3. Physical constraints validation

## References

- ASE Covalent Radii: Cordero et al., Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J

## License

MIT License


## Authors
