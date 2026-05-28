from pathlib import Path
import matplotlib.pyplot as plt
from nmr_setup import setup_nmr_project

# ==========================================
# 1. AUTO-ROUTING & SYSTEM SETUP
# ==========================================
_, PROJECT_NAME, paths = setup_nmr_project(__file__)
import Plot_func as pf

# ==========================================
# 2. I/O HUB
# ==========================================
csv_dir = paths.csv_dir
out_single = paths.out_single
out_overlay = paths.out_overlay
out_grid = paths.out_grid

# %% =======================================
# 3. IMPORT ALL SPECTRA
# ==========================================

base_dir = Path("/User/Documents/NMR_Spectra")

# Data Format (folder, name, contour, color, label.csv(optional))
data = [
    ("250905_Protein_1/700",    "Protein 1",    7e7,    "black",        "Protein_1_pos.csv"),
    ("250905_Protein_2/700",    "Protein 2",    5e8,    "#000000"),
    ("250905_Protein_3/700",    "Protein 3",    6e8,    "tab:blue"),
    ("250905_Protein_4/700",    "Protein 4",    12e8,   [0.839, 0.153, 0.157]),
]

files, file_names, cont, colors, csv_files = pf.parse_plot_data(data, base_dir)
# colors = pf.colormap(files, maps='viridis')  # use a colormap

# %% SELECT YOUR SPECTRA TYPE

# spectra types 15N, 13C, CON, ZOOM, SIZE

spectra_type = "15N"


# %% SETUP PARAMETERS

dpi = 300               # dpi of the figure
save_png = True         # saves svg files (makes it slower)
save_svg = True         # saves svg files (makes it much slower)

xsize = 5               # size of the figure in inch (or of one subplot in the grid-plot)
ysize = 4

lines = 25              # number of contour lines
factor = 1.2            # distance between contour lines

labelpad_x = -8         # distance of the x label to the axis in individual plots
labelpad_y = 9          # distance of the y label to the axis in individual plots

alpha = 0.7             # transparency for overlay
negative = False        # negative conture levels
neg_color = "magenta"   # negative conture level color

first_x = False         # hides the first x-tick label
first_y = False         # hides the first y-tick label
legend = True           # visibility of the legend in overlay plots

# GRIDPLOT PARAMETERS
title_y = 0.88          # distance of the title to the axis in the grid_plot
grid_x = 0.05           # distance of the x label to the axis in the grid_plot
grid_y = 0.06           # distance of the y label to the axis in the grid_plot

# Peaklabel PARAMETERS
peak_labels = True     # toggle visibility of peak labels from CSV files
expand = (1.6, 1.5)     # buffer size around the text to expand
iter_lim = 42           # how many iterations should be used for the label-placement
peak_seed = 42          # random seed - change if you want a different peak_label placement
label_fontsize = 4      # size of the text


# %% SPECTRA TYPE PARAMETERS

if spectra_type == "15N":
    xlim, ylim = (11.2, 6.2), (134, 104)
    x_label, y_label = "$^1$H [ppm]", "$^{15}$N\n[ppm]"
    y_label_grid = "$^{15}$N [ppm]"
    xticks, xminorticks = 1, 0.5
    yticks, yminorticks = 5, 1
    line_width = 0.3

elif spectra_type == "13C":
    xlim, ylim = (1.51, -0.58), (28, 18.2)
    x_label, y_label = "$^1$H [ppm]", "$^{13}$C\n[ppm]"
    y_label_grid = "$^{13}$C [ppm]"
    xticks, xminorticks = 0.5, 0.1
    yticks, yminorticks = 5, 1
    line_width = 0.6

elif spectra_type == "CON":
    xlim, ylim = (178, 168), (140, 108)
    x_label, y_label = "$^{13}$C [ppm]", "$^{15}$N\n[ppm]"
    y_label_grid = "$^{15}$N [ppm]"
    xticks, xminorticks = 2, 0.5
    yticks, yminorticks = 5, 1
    line_width = 0.3

elif spectra_type == "ZOOM":
    xlim, ylim = (8.8, 8), (118, 114)
    x_label, y_label = "$^1$H [ppm]", "$^{15}$N\n[ppm]"
    y_label_grid = "$^{15}$N [ppm]"
    xticks, xminorticks = 0.2, 0.05
    yticks, yminorticks = 2, 0.5
    line_width = 0.6

elif spectra_type == "SIZE":
    xlim, ylim = (8.9, 6.6), (130, 110)
    x_label, y_label = "$^1$H [ppm]", "$^{15}$N\n[ppm]"
    y_label_grid = "$^{15}$N [ppm]"
    xticks, xminorticks = 0.5, 0.25
    yticks, yminorticks = 5, 1
    line_width = 0.3
    xsize = 4
    ysize = 7

else:
    raise ValueError(
        f"Invalid spectra_type: '{spectra_type}'. Must be one of the available types!"
    )


# %% ADDITIONAL SETUP PARAMETERS

# LINEWIDTH AND FONT SIZE
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["svg.fonttype"] = "none"  # set to "path" if you want the text as curves

tick_width = 1
plt.rcParams["xtick.major.width"] = tick_width
plt.rcParams["ytick.major.width"] = tick_width
plt.rcParams["xtick.minor.width"] = tick_width
plt.rcParams["ytick.minor.width"] = tick_width

axis_labels = True
titles = True
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 12
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = 8
plt.rcParams["figure.titlesize"] = 20

# FONT FAMILY
# -- Choose ARIAL or TIMES--
Arial = True
if Arial:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["mathtext.fontset"] = "dejavusans"
else:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["axes.facecolor"] = "None"

# %% LOAD DATA INTO THE MEMORY
print("Loading NMR data into memory...")
dic_all, data_all = pf.read_data(files)

# %% INFO PARAMETERS
# -- Info parameters which will be transferred to the plot functions --
p = {
    "dic_all": dic_all,
    "data_all": data_all,
    "files": files,
    "file_names": file_names,
    "cont": cont,
    "colors": colors,
    "dpi": dpi,
    "xlim": xlim,
    "ylim": ylim,
    "xsize": xsize,
    "ysize": ysize,
    "lines": lines,
    "factor": factor,
    "line_width": line_width,
    "xticks": xticks,
    "xminorticks": xminorticks,
    "yticks": yticks,
    "yminorticks": yminorticks,
    "labelpad_x": labelpad_x,
    "labelpad_y": labelpad_y,
    "title_y": title_y,
    "x_label": x_label,
    "y_label": y_label,
    "y_label_grid": y_label_grid,
    "grid_x": grid_x,
    "grid_y": grid_y,
    "titles": titles,
    "labels": axis_labels,
    "alpha": alpha,
    "negative": negative,
    "neg_color": neg_color,
    "save_svg": save_svg,
    "save_png": save_png,
    "first_x": first_x,
    "first_y": first_y,
    "legend": legend,
    "peak_labels": peak_labels,
    "csv_files": csv_files.tolist(),
    "label_fontsize": label_fontsize,
    "expand": expand,
    "iter_lim": iter_lim,
    "peak_seed": peak_seed,
    "csv_dir": csv_dir,
    "out_single": out_single,
    "out_grid": out_grid,
    "out_overlay": out_overlay,
}


# %% FILES TO OVERLAY
# this gives you the number of the spectrum for the overlay plot

print("\n----- Spectrum Menu -----")
for i, name in enumerate(file_names):
    print(f" [{i}] -> {name}")
print("-------------------------\n")

# %% PLOT ALL SPECTRA
if True:
    pf.plot_everything(p)

# %% PLOT INDIVIDUAL OVERLAYS
# -- you need to define a "name" for the overlay plot --
# -- you can plot all combinations, the last number will be on top --

if True:
    pf.overlay(p, [1, 2], name="over1")
    pf.overlay(p, [2, 1], name="over2")

# %% PLOT ALL SPECTRA IN A GRID
# -- in "grid_plot" you need to define a row and a col for the grid --
if True:
    pf.grid_plot(p, row=2, col=2)

    # %% -- in "grid_plot_over" you need to define a spectra to overlay all others --

    pf.grid_plot_over(p, over=0, row=2, col=2, reverse=False)

    # %% -- in "grid_plot_over_xp" you need to define overlay groups for each subplot --

    pf.grid_plot_over_xp(
        p,
        row=2,
        col=2,
        overlay_groups=[
            [1, 2],
            [2, 1],
            [1, 2, 3],
        ],
    )
