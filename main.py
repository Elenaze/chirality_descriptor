import os
import json
import time
import pandas as pd
from pathlib import Path
import sys
from HCM_atomic_env import HCM_atomic_env_parallel

# Get the directory containing main.py
SCRIPT_DIR = Path(__file__).parent.absolute()


def print_results(system_name, chi_desc, processing_time):
    """
    Print the results of the Chirality Analysis.
    """
    print("\nChirality Descriptor Analysis Results:")
    print("=" * 50)
    print(f"System: {system_name}")
    print(f"Processing time: {processing_time:.2f} seconds")
    print("-" * 30)
    
    for filename, descriptor_value in chi_desc.items():
        print(f"  {filename}: {descriptor_value:.6f}")
    print()


def analyze_single_system(system_path, angle_step, n_processes, alpha, output_dir):
    """
    Analyze chirality descriptors for a single specific system directory.
    """
    system_path = Path(system_path)
    
    # Find .xyz files in the directory
    xyz_files = list(system_path.glob('*.xyz'))
    if not xyz_files:
        raise ValueError(f"No .xyz files found in {system_path}")
    
    # Use the first .xyz file's stem (filename without extension) for naming
    xyz_filename = xyz_files[0].stem
    system_name = xyz_filename  # Use .xyz filename instead of directory name
    
    # Define filename based on .xyz file and parameters
    filename = f"{system_name}_step_{angle_step}_alpha_{alpha}"
    
    # Check if already processed by looking in the CSV file
    csv_file = output_dir / f"{filename}_step{angle_step}.csv"
    if csv_file.exists():
        df_existing = pd.read_csv(csv_file)
        if system_name in df_existing['system'].values:
            print(f"\nSystem '{system_name}' already processed in: {csv_file}")
            
            try:
                response = input("Do you want to reprocess? (y/N): ").strip().lower()
                if response != 'y':
                    print("Skipping reprocessing...")
                    # Show existing results
                    system_results = df_existing[df_existing['system'] == system_name]
                    print(f"\nExisting results for {system_name}:")
                    for _, row in system_results.iterrows():
                        print(f"  {row['filename']}: {row['descriptor_value']:.6f}")
                    return True
                else:
                    # Remove existing entries for this system before reprocessing
                    df_filtered = df_existing[df_existing['system'] != system_name]
                    df_filtered.to_csv(csv_file, index=False)
                    print(f"Removed existing entries for {system_name} from {csv_file}")
            except (EOFError, KeyboardInterrupt):
                print("Skipping reprocessing...")
                return True
    
    start_time = time.time()
    
    try:
        # Run the chirality analysis on this system directory
        chi_desc = HCM_atomic_env_parallel(
            str(system_path), 
            angle_step,
            alpha_value=alpha,
            n_processes=n_processes
        )
        
        end_time = time.time()
        processing_time = end_time - start_time
        
        if chi_desc:
            print(f"Analysis completed")
            
            # Print and save results
            print_results(system_name, chi_desc, processing_time)
            save_results(system_name, chi_desc, output_dir, angle_step, processing_time, alpha, filename)
            
            return True
            
        else:
            raise ValueError(f"No results found for {system_name}")
            
    except Exception as e:
        print(f"Error processing {system_name}: {e}")
        # Save error info to log
        error_log = output_dir / f"errors_{angle_step}.log"
        with open(error_log, 'a') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}: Error in {system_name}: {e}\n")
        return False


def read_config_file(config_path):
    """
    Read configuration parameters from config.txt file.
    Returns: (alpha, angle_step, n_processes, directory_path)
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    params = []
    with open(config_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Remove inline comments
                if '#' in line:
                    line = line.split('#')[0].strip()
                if line:
                    params.append(line)
    
    if len(params) < 4:
        raise ValueError(f"Config file must contain 4 parameters: alpha, angle_step, n_processes, directory_path. Found {len(params)}")
    
    try:
        alpha = float(params[0])
        angle_step = int(params[1])
        n_processes = int(params[2])
        directory_path = params[3]
        
        return alpha, angle_step, n_processes, directory_path
    except (ValueError, IndexError) as e:
        raise ValueError(f"Error parsing config file parameters: {e}")


def save_results(system_name, chi_desc, output_dir, angle_step, processing_time, alpha, filename):
    """Save results for the analyzed system - CSV format only for efficiency and reusability"""
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    # Master CSV file
    csv_data = []
    for file_name, descriptor_value in chi_desc.items():
        csv_data.append({
            'system': system_name,
            'filename': file_name,
            'chirality_value': descriptor_value,
            'alpha': alpha,
            'angle_step': angle_step,
            'processing_time': processing_time,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        })
    
    df_new = pd.DataFrame(csv_data)
    csv_file = output_dir / f"{filename}_step{angle_step}.csv"
    
    # Append to existing CSV or create new one
    if csv_file.exists():
        df_existing = pd.read_csv(csv_file)
        # Filter out any empty or NA columns before concatenation
        df_existing = df_existing.dropna(axis=1, how='all')
        df_new = df_new.dropna(axis=1, how='all')
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
        df_combined.to_csv(csv_file, index=False)
        print(f"Results appended to: {csv_file}")
    else:
        df_new.to_csv(csv_file, index=False)
        print(f"Results saved to: {csv_file}")
    
    # Optional: Also save individual JSON for detailed analysis (if needed)
    # Uncomment the lines below if you want to keep individual system JSONs
    # system_json = output_dir / f"{system_name}_step_{angle_step}.json"
    # system_data = {
    #     'system': system_name,
    #     'angle_step': angle_step,
    #     'alpha': alpha,
    #     'processing_time_seconds': processing_time,
    #     'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    #     'results': chi_desc
    # }
    # with open(system_json, 'w') as f:
    #     json.dump(system_data, f, indent=4)


def main():
    """Main function to run chirality descriptor analysis based on config file"""
    
    # Check for config file
    config_file = SCRIPT_DIR / 'config.txt'
    
    if not config_file.exists():
        print(f"Config file not found: {config_file}")
        print("\nCreate a config.txt file with the following format:")
        print("1.25                     # alpha value for atomic environment analysis")
        print("30                       # angle_step (degrees)")
        print("4                        # n_processes")
        print("data/SR/cAlaAla         # directory with single .xyz file (relative to repository root)")
        sys.exit(1)
    
    try:
        # Read configuration from file
        alpha, angle_step, n_processes, directory_path = read_config_file(config_file)
        
        print(f"Configuration loaded from: {config_file}")
        print(f"Alpha: {alpha}")
        print(f"Angle step: {angle_step}°") 
        print(f"Processes: {n_processes}")
        print(f"Directory: {directory_path}")
        
        # Make directory path absolute if it's relative (relative to script directory)
        input_path = Path(directory_path)
        if not input_path.is_absolute():
            input_path = SCRIPT_DIR / input_path
        
        # Setup output directory
        output_dir = SCRIPT_DIR / 'output'
        output_dir.mkdir(exist_ok=True)
        
        # Check if input path exists and is a directory
        if not input_path.exists():
            raise FileNotFoundError(f"Directory {input_path} does not exist.")
        
        if not input_path.is_dir():
            raise ValueError(f"Path {input_path} is not a directory.")
        
        # Run single system analysis
        chirality = analyze_single_system(input_path, angle_step, n_processes, alpha, output_dir)
         
        return chirality

    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    chirality = main()
    if not chirality:
        print("Analysis failed. Please check the error messages above.")
        sys.exit(1)