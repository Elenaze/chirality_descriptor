from PDF import PDFAnalyzer, AnalysisResults
from datetime import datetime
from pathlib import Path

def main():
    # Get the directory containing main.py
    script_dir = Path(__file__).parent.absolute()
    # Get the root directory (one level up from script_dir)
    root_dir = script_dir.parent
    
    # Create results directory if it doesn't exist
    output_dir = script_dir / "results"
    output_dir.mkdir(exist_ok=True)
    
    # Initialize analyzer - use absolute path to data file
    xyz_file = root_dir / "data" / "combined_minima_all.xyz"
    
    # Check if file exists
    if not xyz_file.exists():
        raise FileNotFoundError(f"Data file not found: {xyz_file}")
    
    print(f"Loading data from: {xyz_file}")
    analyzer = PDFAnalyzer(str(xyz_file))
    
    # Run analysis
    pair_distances = analyzer.analyze()
    
    # Get all results
    radii = analyzer.determine_same_atom_radii()
    radii = analyzer.determine_different_atom_radii(radii)
    X_max, X_min = analyzer.find_max_min(min_peak_height=0.002)
    
    # Clean up radii
    for atom, radius_list in radii.items():
        radii[atom] = next((r for r in radius_list if r != 0), 0)
    
    # Get alpha values
    df_alpha = analyzer.determine_alpha(radii)
    
    # Create summary data
    pair_summary = {}
    for (atom1, atom2), distances in pair_distances.items():
        pair_summary[(atom1, atom2)] = {
            'count': len(distances),
            'min_distance': min(distances),
            'max_distance': max(distances)
        }
    
    # Create results object
    results = AnalysisResults(
        timestamp=datetime.now().isoformat(),
        xyz_file=xyz_file,
        pair_distances=pair_distances,
        kde_data=analyzer.kde_data,
        X_max=X_max,
        X_min=X_min,
        radii=radii,
        alpha_values=df_alpha,
        pair_summary=pair_summary
    )
    
    # Save results
    results.save(output_dir)
    
    # Print summary
    print("\nAnalysis Results Summary")
    print("=" * 40)
    print(f"\nPair Distance Statistics:")
    for (atom1, atom2), stats in pair_summary.items():
        print(f"\n{atom1}-{atom2} pairs:")
        print(f"  Count: {stats['count']}")
        print(f"  Range: {stats['min_distance']:.3f} - {stats['max_distance']:.3f} Å")
    
    print("\nAtomic Radii:")
    for atom, radius in radii.items():
        print(f"  {atom}: {radius:.3f} Å")
    
    print("\nAlpha Values:")
    print(df_alpha)
    
    print(f"\nResults saved to: {output_dir}")
     # Create plots subdirectory
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    analyzer.plot_pair_distributions(plots_dir)

if __name__ == "__main__":
    main()