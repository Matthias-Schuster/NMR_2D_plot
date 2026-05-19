# NMR_2D_plot 📈

A streamlined, highly customizable Python pipeline for visualizing and rendering publication-quality 2D NMR spectra. 

Built on `matplotlib` and `nmrglue`, this tool reads Bruker processed data (`pdata`) directly and generates single spectra, overlays, and complex grid plots.

## ✨ Features
* **Direct Bruker Integration:** Reads processed NMR data natively using `nmrglue`.
* **Smart Contour Caching:** Caches calculated contour vertices in memory, drastically speeding up the rendering of complex overlays and grid plots.
* **Versatile Plot Types:** Generate individual plots, multi-spectrum overlays, multi-panel grids, and overlaid grids.
* **Pre-configured Spectra Presets:** Out-of-the-box styling and axis limits for standard 2D experiments (e.g., 15N-HSQC, 13C-HSQC, CON).
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
```

## 🚀 Quickstart & Usage
The main user interface is the NMR_plot.py script. You do not need to interact with Plot_func.py unless you want to change the underlying matplotlib logic.

### 1. Define Your Data
In NMR_plot.py, locate the data list. Add your spectra using the following tuple format: ('folder_path', 'Spectrum_Name', contour_base_level, 'color', 'optional_label.csv').
The nef_extract script can generate a template for this input if you have loaded the spectra into ccpnmr3.

```Python
base_dir = Path('/User/Documents/NMR_Spectra/')

data = [
    # (Folder, Name, Contour Start, Color, CSV File)
    ('250905_Protein_1/700', "Protein 1", 7e7, "black", "Protein1_labels.csv"),
    ('250905_Protein_2/700', "Protein 2", 5e8, "tab:blue"), # CSV is optional!
]
```

### 2. Choose Your Spectrum Type
Select the pre-configured axis limits and labels for your experiment:

```Python
spectra_type = "15N"  # Options: "15N", "13C", "CON", "ZOOM", "SIZE"
```

### 3. Generate Plots
At the bottom of NMR_plot.py, toggle the boolean flags (if True:) to generate the plot types you need. The script will automatically create the necessary output folders (e.g., results/single, results/over).
* Single Plots: Plots each spectrum individually.
* Overlays: Specify indices from your data list to overlay (e.g., [0, 1] overlays Protein 1 and Protein 2). The last index is plotted on top.
* Grid Plots: Plots all loaded spectra in an $M \times N$ grid.
* Overlaid Grids: Select a reference spectrum to overlay across all panels in your grid.
* Custom Grid Overlays: Define specific groupings of spectra for each individual grid panel using grid_plot_over_xp.

## 🏷️ Peak Labels via CSV (Optional)
You can automatically annotate your NMR spectra with peak labels by providing a CSV file for each spectrum.

### 1. Directory Setup
Place your label CSV files inside a folder named csv/ in the same directory where you run your plotting script.

```plaintext
├── NMR_plot.py
├── Plot_func.py
└── csv/
    ├── Protein1_labels.csv
    └── Protein2_labels.csv
```

### 2. CSV Format
The CSV file must contain the following exactly named columns:

* Residue: The text label to display (e.g., "A24", "ASF1A").
* position_1: The X-axis coordinate (e.g., ¹H ppm).
* position_2: The Y-axis coordinate (e.g., ¹⁵N or ¹³C ppm).

Example Protein1_labels.csv:

```plaintext
Residue,position_1,position_2
A24,8.15,123.4
G50,8.42,110.2
```

Once your CSV is in the csv/ folder, just add the filename as the 5th element in your data tuple (as shown in Step 1) and the script will automatically render and space out your labels!

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
