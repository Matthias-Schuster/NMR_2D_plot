# NMR_2D_plot 📈

A streamlined, highly customizable Python pipeline for visualizing and rendering publication-quality 2D NMR spectra. 

Built on `matplotlib` and `nmrglue`, this tool reads Bruker processed data (`pdata`) directly and generates single spectra, overlays, and complex grid plots.

## ✨ Features
* **Direct Bruker Integration:** Reads processed NMR data natively using `nmrglue`.
* **Smart Contour Caching:** Caches calculated contour vertices in memory, drastically speeding up the rendering of complex overlays and grid plots.
* **Versatile Plot Types:** Generate individual plots, multi-spectrum overlays, multi-panel grids, and overlaid grids.
* **Automated Peak Labeling:** Import peak positions directly from a CSV file with an intelligent layout algorithm (`adjustText`) to prevent label overlapping.
* **Auto-Routing Architecture:** Effortlessly manage multiple projects. The pipeline automatically routes your inputs and outputs into dedicated project folders based on your script's name.
* **Publication Ready:** Easy export to high-DPI PNGs or vector graphic SVGs.

## ⚙️ Installation
We recommend using `conda` to manage dependencies. An `environment.yml` file is provided for quick setup.

```bash
# Clone the repository
git clone https://github.com/Matthias-Schuster/NMR_2D_plot.git
cd NMR_2D_plot

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate NMR_2D_plot


# Alternatively you can compare or update your environment 
conda compare -n your_env_name environment.yml
conda env update -n your_env_name -f environment.yml
```

## 📁 Project Architecture
This pipeline uses a strict Input/Output hub to keep your data organized. The core engine (Plot_func.py) stays at the root, while your individual plotting scripts live in plot_scripts/.

```plaintext
NMR_2D_plot/
├── Plot_func.py                 # Core plotting engine!
├── environment.yml
│
├── plot_scripts/                # ⬅️ Work from here!
│   ├── nmr_setup.py             # Auto-routing script
│   └── NMR_plot_Template.py     # Master template
│
├── input_csv/                   # (Auto-generated) Peaklist inputs
└── results/                     # (Auto-generated) Plot outputs
```


## 🚀 Quickstart & Usage

### 1. Create a New Project Script
Navigate to the plot_scripts/ folder. Duplicate NMR_plot_Template.py and rename it for your specific project. For example: NMR_plot_MyProtein.py.

The pipeline will automatically read the name MyProtein and create dedicated input_csv/MyProtein/ and results/MyProtein/ folders for you the first time you run it!

### 2. Define Your Data
Open your newly renamed script. In the data block, add your spectra using the following tuple format: ('folder_path', 'Spectrum_Name', contour_base_level, 'color', 'optional_label.csv').
The nef_extract script can generate a template for this input if you have loaded the spectra into ccpnmr3.

```Python
base_dir = Path('/User/Documents/NMR_Spectra/')

data = [
    # (Folder, Name, Contour Start, Color, CSV File)
    ('250905_Protein_1/700', "Protein 1", 7e7, "black", "Protein1_labels.csv"),
    ('250905_Protein_2/700', "Protein 2", 5e8, "tab:blue"), # CSV is optional!
    ('250905_Protein_3/700', "Protein 3", 6e8, ["tab:orange", "tab:cyan"]) # a second color activates negative contours!
]
```

### 3. Choose Your Spectrum Type
Select the pre-configured axis limits and labels for your experiment:

```Python
pf.set_style(p, "15N")  # Options: "15N", "13C", "CON", "ZOOM", "IDP"
```
You can modify the Styles of the SPECTRA_STYLES dictionary!

### 4. Generate Plots
At the bottom of your script, toggle the boolean flags (if True:) to generate the plot types you need. The script will automatically save them to your project's results folder.
* Single Plots: Plots each spectrum individually.
* Overlays: Specify indices from your data list to overlay (e.g., [0, 1] overlays Protein 1 and Protein 2). The last index is plotted on top.
* Grid Plots: Plots all loaded spectra in an $M \times N$ grid.
* Overlaid Grids: Select a reference spectrum to overlay across all panels in your grid.
* Custom Grid Overlays: Define specific groupings of spectra for each individual grid panel using grid_plot_over_xp.

For example:

```Python
pf.plot_everything(p)

pf.overlay(p, [1], name="prot1")
pf.overlay(p, [1, 2], name="over1")
pf.overlay(p, [2, 1], name="over2")

pf.grid_plot(p, row=2, col=2)
```

### 5. Generate Slices from 3D spectra

* plot_strip: Plots a single slice of a 3D spectrum (or overlays of 3D spectra).
* strip_grid_plot: Plots multiple slices of a 3D spectrum (or overlays of 3D spectra).

Slices can be defined like this

```Python
my_slices = [
    {"z_slice": 120.75, "x_slice": 8.23, "title": "Asn12"},
    {"z_slice": 121.38, "x_slice": 8.08, "title": "Glu13"},
    {"z_slice": 122.47, "x_slice": 7.85, "title": "Val14"},
]

pf.strip_grid_plot(
    p,
    [4, 5],
    slices=my_slices,
    slice_axis=1,
    row=1,
    col=3,
    strip_width=0.2,
    name="Sequential_Walk",
    reverse=False,
)
```

## 🏷️ Peak Labels via CSV (Optional)
You can automatically annotate your NMR spectra with peak labels by providing a CSV file for each spectrum (only for 2D spectra).

### 1. Directory Setup
Once you run your project script (e.g., NMR_plot_MyProtein.py) for the first time, it will generate an empty folder for your labels at input_csv/MyProtein/. Place your CSV files in there.

### 2. CSV Format
The CSV file must contain the following exactly named columns:

* Residue: The text label to display (e.g., "A24").
* position_1: The X-axis coordinate (e.g., ¹H ppm).
* position_2: The Y-axis coordinate (e.g., ¹⁵N or ¹³C ppm).

Example Protein1_labels.csv:

```plaintext
Residue,position_1,position_2
A24,8.15,123.4
G50,8.42,110.2
```

Once your CSV is in the input_csv/MyProtein/ folder, just ensure the filename matches the 5th element in your data tuple (as shown in Step 2) and the script will automatically render and space out your labels!

## 🎨 Customization
You can fine-tune your figures in the SETUP PARAMETERS section of NMR_plot.py:
* Contour adjustments: Change lines (number of contour levels) and factor (spacing between levels).
* Negative contours: Set negative = True and pick a neg_color to render negative intensities as dashed lines.
* Typography: Toggle the Arial flag to easily switch between sans-serif (Arial) and serif (Times New Roman) fonts.
* Exporting: Set save_png = True or save_svg = True. (Note: SVG generation is significantly slower for contour-heavy NMR spectra).

## 🤝 Contributing
Contributions, issues, and feature requests are welcome!

## 📝 License
Distributed under the MIT License.
