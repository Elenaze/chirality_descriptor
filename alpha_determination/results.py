from dataclasses import dataclass
from datetime import datetime
import json
import pickle
from pathlib import Path
import pandas as pd
from typing import Dict, Tuple, List

@dataclass
class AnalysisResults:
    """Container for PDF analysis results"""
    timestamp: str
    xyz_file: str
    pair_distances: Dict[Tuple[str, str], List[float]]
    kde_data: Dict[Tuple[str, str], pd.DataFrame]
    X_max: Dict[Tuple[str, str], float]
    X_min: Dict[Tuple[str, str], float]
    radii: Dict[str, float]
    alpha_values: pd.DataFrame
    pair_summary: Dict[Tuple[str, str], Dict[str, float]]
    
    def save(self, output_dir: str):
        """Save results to directory"""
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        
        # Convert alpha values DataFrame to dict
        alpha_dict = {
            f"{row['Pair'][0]}-{row['Pair'][1]}": row['Alpha'] 
            for _, row in self.alpha_values.iterrows()
        }
        
        # Create comprehensive analysis results
        analysis_results = {
            'metadata': {
                'timestamp': self.timestamp,
                'xyz_file': str(self.xyz_file)
            },
            'atomic_radii': self.radii,
            'maxima': {f"{k[0]}-{k[1]}": v for k, v in self.X_max.items()},
            'minima': {f"{k[0]}-{k[1]}": v for k, v in self.X_min.items()},
            'pair_summary': {f"{k[0]}-{k[1]}": v for k, v in self.pair_summary.items()},
            'alpha_values': alpha_dict
        }
        
        # Save complete analysis results
        with open(out_path / 'analysis_results.json', 'w') as f:
            json.dump(analysis_results, f, indent=2)
        
        # Save KDE data separately as they're large DataFrames
        for key, df in self.kde_data.items():
            df.to_csv(out_path / f'kde_{key[0]}_{key[1]}.csv')
        
        # Save full results object for programmatic access
        with open(out_path / 'full_results.pkl', 'wb') as f:
            pickle.dump(self, f)