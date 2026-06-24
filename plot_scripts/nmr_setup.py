from pathlib import Path
import sys


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

    # Make root directory modules importable
    root_dir_str = str(root_dir)
    if root_dir_str not in sys.path:
        sys.path.insert(0, root_dir_str)

    # Extract project name from filename
    script_name = script_path.stem

    if script_name.startswith(prefix):
        project_name = script_name[len(prefix):]
    else:
        project_name = script_name

    # Define paths
    project_out = root_dir / results_folder / project_name

    # Group paths in a dictionary
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

    return project_name, paths


def build_plot_dict(**kwargs):
    """
    Builds the 'p' dictionary from the local variables of the main script.
    """
    paths = kwargs.get("paths", {})

    p = {
        "spectra_type": kwargs.get("spectra_type"),
        "dic_all": kwargs.get("dic_all"),
        "data_all": kwargs.get("data_all"),
        "files": kwargs.get("files"),
        "file_names": kwargs.get("file_names"),
        "cont": kwargs.get("cont"),
        "colors": kwargs.get("colors"),
        # Pull parameters determined by your spectra_type if/elif block
        "xlim": kwargs.get("xlim"),
        "ylim": kwargs.get("ylim"),
        "x_label": kwargs.get("x_label"),
        "y_label": kwargs.get("y_label"),
        "y_label_grid": kwargs.get("y_label_grid"),
        "xticks": kwargs.get("xticks"),
        "xminorticks": kwargs.get("xminorticks"),
        "yticks": kwargs.get("yticks"),
        "yminorticks": kwargs.get("yminorticks"),
        "line_width": kwargs.get("line_width"),
        # Pull all your other manual settings
        "dpi": kwargs.get("dpi", 300),
        "xsize": kwargs.get("xsize", 5),
        "ysize": kwargs.get("ysize", 4),
        "lines": kwargs.get("lines", 25),
        "factor": kwargs.get("factor", 1.2),
        "labelpad_x": kwargs.get("labelpad_x", -8),
        "labelpad_y": kwargs.get("labelpad_y", 9),
        "title_y": kwargs.get("title_y", 0.88),
        "grid_x": kwargs.get("grid_x", 0.05),
        "grid_y": kwargs.get("grid_y", 0.06),
        "titles": kwargs.get("titles", True),
        "labels": kwargs.get("axis_labels", True),
        "alpha": kwargs.get("alpha", 0.7),
        "negative": kwargs.get("negative", False),
        "neg_color": kwargs.get("neg_color", "magenta"),
        "save_svg": kwargs.get("save_svg", False),
        "save_png": kwargs.get("save_png", True),
        "first_x": kwargs.get("first_x", True),
        "first_y": kwargs.get("first_y", True),
        "legend": kwargs.get("legend", True),
        "peak_labels": kwargs.get("peak_labels", True),
        "csv_files": kwargs.get("csv_files", []),
        "label_fontsize": kwargs.get("label_fontsize", 4),
        "csv_dir": paths.get("csv_dir", ""),
        "out_single": paths.get("out_single", ""),
        "out_grid": paths.get("out_grid", ""),
        "out_overlay": paths.get("out_overlay", ""),
        "peak_ellipse_x": kwargs.get("peak_ellipse_x", 0.08),
        "peak_ellipse_y": kwargs.get("peak_ellipse_y", 0.8),
        "show_ellipse": kwargs.get("show_ellipse", True),
        "starting_angle": kwargs.get("starting_angle", 0.0),
        "random_seed": kwargs.get("random_seed", 42),
        "cycles": kwargs.get("cycles", 10),
        "iterations": kwargs.get("iterations", 10),
        "connector_width": kwargs.get("connector_width", 0.5),
    }

    return p
