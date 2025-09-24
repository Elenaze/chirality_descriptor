import os
import json
import time
import pandas as pd
from pathlib import Path
import sys
from HCM_atomic_env import HCM_atomic_env_parallel

# Get the directory containing main.py
SCRIPT_DIR = Path(__file__).parent.absolute()
# Get the root directory (one level up from SCRIPT_DIR)
ROOT_DIR = SCRIPT_DIR.parent

def get_completed_systems(output_dir, angle_step):
    """Check which systems have already been completed"""
    completed = {'SR': set(), 'SS': set()}
    
    # Check individual JSON files
    for json_file in output_dir.glob(f"*_*_step_{angle_step}.json"):
        parts = json_file.stem.split('_')
        if len(parts) >= 3:
            chirality = parts[0]
            system_name = '_'.join(parts[1:-2])  # Handle system names with underscores
            if chirality in completed:
                completed[chirality].add(system_name)
                
    return completed

def print_summary(results):
    """
    Print a summary of the chirality analysis results.
    """
    print("Chirality Descriptor Analysis Summary:")
    print("=" * 50)
    
    for chirality, systems in results.items():
        print(f"\nChirality: {chirality}")
        print("-" * 30)
        
        for system_name, files in systems.items():
            print(f"  System: {system_name}")
            for filename, descriptor_value in files.items():
                print(f"    {filename}: {descriptor_value:.6f}")
        print()

def update_progress_file(output_dir, angle_step, chirality, system_name, system_time):
    """Update the progress tracking file"""
    progress_file = output_dir / f"progress_step_{angle_step}.txt"
    
    with open(progress_file, 'a') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}: Completed {chirality}/{system_name} in {system_time:.2f}s\n")

def save_system_results(system_name, chirality, chi_desc, output_dir, angle_step, system_time):
    """Save results for a single system immediately after processing"""
    
    # Individual system JSON file
    system_json = output_dir / f"{chirality}_{system_name}_step_{angle_step}.json"
    system_data = {
        'chirality': chirality,
        'system': system_name,
        'angle_step': angle_step,
        'processing_time_seconds': system_time,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'results': chi_desc
    }
    
    with open(system_json, 'w') as f:
        json.dump(system_data, f, indent=4)
    print(f"    Individual results saved to: {system_json}")
    
    # Append to chirality-specific DAT file
    dat_file = output_dir / f"{angle_step}_step_env_hd_{chirality}.dat"
    with open(dat_file, 'a') as f:  # Append mode
        for filename, descriptor_value in chi_desc.items():
            f.write(f"{filename}\t{descriptor_value:.6f}\t{system_name}\n")
    
    # Append to master CSV file
    csv_file = output_dir / f"chirality_descriptors_step_{angle_step}.csv"
    csv_data = []
    for filename, descriptor_value in chi_desc.items():
        csv_data.append({
            'chirality': chirality,
            'system': system_name,
            'filename': filename,
            'descriptor_value': descriptor_value,
            'processing_time': system_time,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        })
    
    df_new = pd.DataFrame(csv_data)
    
    # Append to existing CSV or create new one
    if csv_file.exists():
        df_existing = pd.read_csv(csv_file)
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
        df_combined.to_csv(csv_file, index=False)
    else:
        df_new.to_csv(csv_file, index=False)

def save_final_summary(results, output_dir, angle_step, total_runtime):
    """Save final comprehensive summary after all processing is complete"""
    
    # Final comprehensive JSON
    json_file = output_dir / f"chirality_descriptors_complete_step_{angle_step}.json"
    export_dict = {
        'parameters': {
            'angle_step': angle_step,
            'total_runtime_seconds': total_runtime,
            'completion_timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        },
        'results': results
    }
    
    with open(json_file, 'w') as f:
        json.dump(export_dict, f, indent=4)
    print(f"\nFinal comprehensive results saved to: {json_file}")
    
    # Add timing summary to DAT files
    for chirality in ['SR', 'SS']:
        dat_file = output_dir / f"{angle_step}_step_env_hd_{chirality}.dat"
        if dat_file.exists():
            with open(dat_file, 'a') as f:
                f.write(f"\n# Total_runtime_seconds\t{total_runtime:.4f}\n")
                f.write(f"# Completion_timestamp\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")

def analyze_single_system(system_path, angle_step, n_processes, output_dir):
    """
    Analyze chirality descriptors for a single specific system directory.
    """
    system_path = Path(system_path)
    system_name = system_path.name
    
    # Determine chirality from parent directory name or path
    chirality = None
    if 'SR' in str(system_path):
        chirality = 'SR'
    elif 'SS' in str(system_path):
        chirality = 'SS'
    else:
        # Try to infer from parent directory
        parent_name = system_path.parent.name
        if parent_name in ['SR', 'SS']:
            chirality = parent_name
        else:
            raise ValueError(f"Cannot determine chirality from path: {system_path}")
    
    print(f"Single System Chirality Descriptor Analysis")
    print(f"==========================================")
    print(f"System path: {system_path}")
    print(f"System name: {system_name}")
    print(f"Chirality: {chirality}")
    print(f"Angle step: {angle_step}°")
    print(f"Processes: {n_processes}")
    
    # Check if already completed
    system_json = output_dir / f"{chirality}_{system_name}_step_{angle_step}.json"
    if system_json.exists():
        print(f"\n⚠️  System already processed! Found existing results at:")
        print(f"   {system_json}")
        
        # Check if running interactively or from config file
        try:
            if sys.stdin.isatty():  # Interactive terminal
                response = input("Do you want to reprocess? (y/N): ").strip().lower()
                if response != 'y':
                    print("Skipping reprocessing. Loading existing results...")
                    with open(system_json, 'r') as f:
                        system_data = json.load(f)
                    return {chirality: {system_name: system_data['results']}}, system_data['processing_time_seconds']
            else:  # Non-interactive (config file or pipe)
                print("Non-interactive mode: skipping reprocessing. Loading existing results...")
                with open(system_json, 'r') as f:
                    system_data = json.load(f)
                return {chirality: {system_name: system_data['results']}}, system_data['processing_time_seconds']
        except:
            # Fallback: skip reprocessing
            print("Skipping reprocessing. Loading existing results...")
            with open(system_json, 'r') as f:
                system_data = json.load(f)
            return {chirality: {system_name: system_data['results']}}, system_data['processing_time_seconds']
    
    # Initialize output files if needed
    dat_file = output_dir / f"{angle_step}_step_env_hd_{chirality}.dat"
    if not dat_file.exists():
        with open(dat_file, 'w') as f:
            f.write(f"# Chirality Descriptor Results - {chirality}\n")
            f.write(f"# Angle step: {angle_step} degrees\n")
            f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Format: filename\tdescriptor_value\tsystem_name\n")
            f.write("#" + "="*50 + "\n")
    
    print(f"\n🚀 Starting analysis...")
    start_time = time.time()
    
    try:
        # Run the chirality analysis on this system directory
        chi_desc = HCM_atomic_env_parallel(
            str(system_path), 
            angle_step, 
            n_processes=n_processes
        )
        
        end_time = time.time()
        system_time = end_time - start_time
        
        if chi_desc:
            print(f"✅ Analysis completed successfully!")
            print(f"   Processing time: {system_time:.2f} seconds")
            print(f"   Found {len(chi_desc)} structures")
            
            # Save results immediately
            save_system_results(system_name, chirality, chi_desc, output_dir, angle_step, system_time)
            update_progress_file(output_dir, angle_step, chirality, system_name, system_time)
            
            results = {chirality: {system_name: chi_desc}}
            return results, system_time
            
        else:
            raise ValueError(f"No results found for {system_name}")
            
    except Exception as e:
        print(f"❌ Error processing {system_name}: {e}")
        # Save error info to log
        error_log = output_dir / f"errors_{angle_step}.log"
        with open(error_log, 'a') as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}: Error in {chirality}/{system_name}: {e}\n")
        raise

def analyze_chirality_data(data_path, angle_step, n_processes, output_dir):
    """
    Analyze chirality descriptors for all systems in the data directory.
    Saves results incrementally after each system is processed.
    """
    results = {}
    total_start_time = time.time()
    
    # Check for existing results and potentially resume
    completed_systems = get_completed_systems(output_dir, angle_step)
    total_completed = sum(len(systems) for systems in completed_systems.values())
    
    if total_completed > 0:
        print(f"\n🔄 Found {total_completed} already completed systems:")
        for chirality, systems in completed_systems.items():
            if systems:
                print(f"  {chirality}: {', '.join(sorted(systems))}")
        print("Will skip these systems and continue with remaining ones.\n")
    
    # Initialize output files (append mode if resuming)
    for chirality in ['SR', 'SS']:
        dat_file = output_dir / f"{angle_step}_step_env_hd_{chirality}.dat"
        if not dat_file.exists():
            with open(dat_file, 'w') as f:
                f.write(f"# Chirality Descriptor Results - {chirality}\n")
                f.write(f"# Angle step: {angle_step} degrees\n")
                f.write(f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"# Format: filename\tdescriptor_value\tsystem_name\n")
                f.write("#" + "="*50 + "\n")
    
    for chirality in ['SR', 'SS']:
        chirality_path = data_path / chirality
        print(f"\nProcessing chirality: {chirality}")
        print(f"Looking in folder: {chirality_path}")
        
        if not chirality_path.exists():
            print(f"Folder does not exist: {chirality_path}")
            continue
            
        results[chirality] = {}
        
        # Find all system directories
        system_dirs = [d for d in chirality_path.iterdir() if d.is_dir()]
        print(f"Found {len(system_dirs)} system directories")
        
        for system_dir in system_dirs:
            system_name = system_dir.name
            
            # Skip if already completed
            if system_name in completed_systems.get(chirality, set()):
                print(f"\n  ⏭️  Skipping {system_name} (already completed)")
                # Load existing results for final summary
                json_file = output_dir / f"{chirality}_{system_name}_step_{angle_step}.json"
                if json_file.exists():
                    with open(json_file, 'r') as f:
                        system_data = json.load(f)
                        results[chirality][system_name] = system_data['results']
                continue
            
            print(f"\n  Processing system: {system_name}")
            start_time = time.time()
            
            try:
                # Run the chirality analysis on this system directory
                chi_desc = HCM_atomic_env_parallel(
                    str(system_dir), 
                    angle_step, 
                    n_processes=n_processes
                )
                
                end_time = time.time()
                system_time = end_time - start_time
                
                if chi_desc:
                    results[chirality][system_name] = chi_desc
                    print(f"    Completed in {system_time:.2f} seconds")
                    print(f"    Found {len(chi_desc)} structures")
                    
                    # ✅ SAVE RESULTS IMMEDIATELY
                    save_system_results(system_name, chirality, chi_desc, output_dir, angle_step, system_time)
                    
                    # Update progress tracking
                    update_progress_file(output_dir, angle_step, chirality, system_name, system_time)
                    
                else:
                    print(f"    No results found for {system_name}")
                    
            except Exception as e:
                print(f"    Error processing {system_name}: {e}")
                # Save error info to log
                error_log = output_dir / f"errors_{angle_step}.log"
                with open(error_log, 'a') as f:
                    f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}: Error in {chirality}/{system_name}: {e}\n")
                continue
    
    total_end_time = time.time()
    total_runtime = total_end_time - total_start_time
    
    return results, total_runtime

def main():
    """Main function to run chirality descriptor analysis"""
    
    # Check if we have command line arguments or should read from stdin
    if len(sys.argv) >= 2:
        # Parse command line arguments
        angle_step = int(sys.argv[1])
        
        if len(sys.argv) >= 3:
            input_path = Path(sys.argv[2])
            if not input_path.is_absolute():
                input_path = ROOT_DIR / input_path
        else:
            input_path = ROOT_DIR / 'data'
        
        n_processes = int(sys.argv[3]) if len(sys.argv) >= 4 else 8
        
        # Check for single system mode
        if len(sys.argv) >= 5:
            system_directory = Path(sys.argv[4])
            if not system_directory.is_absolute():
                system_directory = ROOT_DIR / system_directory
            single_system_mode = True
        else:
            system_directory = None
            single_system_mode = False
            data_path = input_path
        
    else:
        # Try to read from stdin (for input redirection)
        try:
            print("Reading parameters from stdin...")
            lines = []
            for line in sys.stdin:
                line = line.strip()
                if line and not line.startswith('#'):  # Skip empty lines and comment-only lines
                    # Remove inline comments
                    if '#' in line:
                        line = line.split('#')[0].strip()
                    if line:  # Check again after removing comments
                        lines.append(line)
            
            if len(lines) < 1:
                raise ValueError("No parameters provided")
            
            angle_step = int(lines[0])
            
            if len(lines) >= 2:
                input_path = Path(lines[1])
                if not input_path.is_absolute():
                    input_path = ROOT_DIR / input_path
            else:
                input_path = ROOT_DIR / 'data'
            
            n_processes = int(lines[2]) if len(lines) >= 3 else 8
            
            # Check if input_path is a specific system directory or general data directory
            if len(lines) >= 4:
                system_directory = Path(lines[3])
                if not system_directory.is_absolute():
                    system_directory = ROOT_DIR / system_directory
                single_system_mode = True
            else:
                system_directory = None
                single_system_mode = False
                data_path = input_path
            
        except (ValueError, IndexError) as e:
            print("Usage: python main.py <angle_step> [data_directory] [n_processes] [system_directory]")
            print("Or use input redirection: python main.py < config.txt")
            print("\nConfig file format:")
            print("  Line 1: angle_step (required)")
            print("  Line 2: data_directory (optional, default: ../data)")
            print("  Line 3: n_processes (optional, default: 8)")
            print("  Line 4: system_directory (optional, for single system analysis)")
            print("\nExamples:")
            print("  python main.py 30")
            print("  python main.py 30 ../data 8")
            print("  python main.py 30 ../data 8 ../data/SR/cAlaAla")
            print(f"Error reading parameters: {e}")
            sys.exit(1)
    
    # Setup output directory
    output_dir = SCRIPT_DIR / 'output'
    output_dir.mkdir(exist_ok=True)
    
    try:
        if single_system_mode:
            # Single system analysis mode
            if not system_directory.is_dir():
                raise FileNotFoundError(f"System directory {system_directory} does not exist.")
            
            # Run single system analysis
            results, total_runtime = analyze_single_system(system_directory, angle_step, n_processes, output_dir)
            
        else:
            # Multiple systems analysis mode (original behavior)
            print(f"Chirality Descriptor Analysis - Multiple Systems")
            print(f"===============================================")
            print(f"Data directory: {data_path}")
            print(f"Angle step: {angle_step}°")
            print(f"Processes: {n_processes}")
            print(f"Output directory: {output_dir}")
            
            # Validate input directory
            if not data_path.is_dir():
                raise FileNotFoundError(f"Directory {data_path} does not exist.")
            
            # Run analysis with incremental saving
            results, total_runtime = analyze_chirality_data(data_path, angle_step, n_processes, output_dir)
        
        if not results or all(not systems for systems in results.values()):
            raise ValueError("No results found in the analysis.")
        
        # Process and save final summary (for both modes)
        print_summary(results)
        if not single_system_mode:  # Only save comprehensive summary for multiple systems
            save_final_summary(results, output_dir, angle_step, total_runtime)
        
        print(f"\nAnalysis completed successfully!")
        print(f"Total runtime: {total_runtime:.2f} seconds")
        
        return True

    except FileNotFoundError as e:
        print(f"File system error: {str(e)}")
        return False
    except ValueError as e:
        print(f"Value error: {str(e)}")
        return False
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        return False

if __name__ == "__main__":
    success = main()
    if not success:
        print("Analysis failed. Please check the error messages above.")
        sys.exit(1)