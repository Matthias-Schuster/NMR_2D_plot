from pathlib import Path
import sys
import importlib


def setup_nmr_project(
    script_file,
    prefix="NMR_plot_",
    results_folder="results",
    csv_folder="input_csv",
    create_output_dirs=True,
    verbose=True,
):
    """
    Automatic setup for NMR plotting scripts.

    It determines:
    - where the script is located
    - where the project root is
    - the project name from the script filename
    - where results should be saved
    - where project-specific CSV files are stored
    """

    # Location of the currently running plot script
    script_path = Path(script_file).resolve()
    script_dir = script_path.parent

    # Project root is one level above plot_scripts/
    root_dir = script_dir.parent

    # Make Plot_func.py importable from the project root
    root_dir_str = str(root_dir)
    if root_dir_str not in sys.path:
        sys.path.insert(0, root_dir_str)

    # Import Plot_func.py from the root directory
    pf = importlib.import_module("Plot_func")

    # Extract project name from filename
    script_name = script_path.stem

    if script_name.startswith(prefix):
        project_name = script_name[len(prefix) :]
    else:
        project_name = script_name

    # Define paths
    project_out = root_dir / results_folder / project_name

    # Group paths in a dictionary instead of a class
    paths = {
        "script_dir": script_dir,
        "root_dir": root_dir,
        "project_out": project_out,
        "out_single": project_out / "single",
        "out_overlay": project_out / "overlay",
        "out_grid": project_out / "grid",
        "csv_dir": root_dir / csv_folder / project_name,
    }

    # Create output folders
    if create_output_dirs:
        paths["out_single"].mkdir(parents=True, exist_ok=True)
        paths["out_overlay"].mkdir(parents=True, exist_ok=True)
        paths["out_grid"].mkdir(parents=True, exist_ok=True)
        paths["csv_dir"].mkdir(parents=True, exist_ok=True)

    if verbose:
        print("\n----- NMR Project Setup -----")
        print(f"Project name: {project_name}")
        print(f"Script dir:   {paths['script_dir']}")
        print(f"Root dir:     {paths['root_dir']}")
        print(f"CSV dir:      {paths['csv_dir']}")
        print(f"Output dir:   {paths['project_out']}")
        print("-----------------------------\n")

    return pf, project_name, paths
