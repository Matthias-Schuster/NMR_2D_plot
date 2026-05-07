from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import Plot_func as pf


# %% IMPORT ALL SPECTRA

# directory of your NMR files
base_dir = Path("/User/Documents/NMR_Spectra/")

# Data Format (folder, name, contour, color)
data = [
    ("250905_Protein_1/700",    "Protein 1",    7e7,    "black"),
    ("250905_Protein_2/700",    "Protein 2",    5e8,    "#000000"),
    ("250905_Protein_3/700",    "Protein 3",    6e8,    "tab:blue"),
    ("250905_Protein_4/700",    "Protein 4",    12e8,   [0.839, 0.153, 0.157]),
]

# subdirectory of your processed files
sub_dir = "pdata/1"

df = pd.DataFrame(data, columns=["folder", "name", "contour", "color"])
df["full_path"] = df["folder"].apply(lambda x: str(base_dir / x / sub_dir))

files = df["full_path"]
file_names = df["name"]
cont = df["contour"]
colors = df["color"]
# colors = pf.colormap(files, maps='viridis')  # use a colormap

# %% SELECT YOUR SPECTRA TYPE

# spectra types 15N, 13C, CON, ZOOM, SIZE

spectra_type = "15N"


# %% SETUP PARAMETERS

dpi = 300               # dpi of the figure
save_png = False        # saves png files (makes it slower)
save_svg = False        # saves svg files (makes it much slower)

xsize = 5               # size of the figure in inch (or of one subplot in the grid-plot)
ysize = 4

lines = 25              # number of contour lines
factor = 1.2            # distance between contour lines

labelpad_x = -8         # distance of the x label to the axis in individual plots
labelpad_y = 9          # distance of the y label to the axis in individual plots

alpha = 0.8             # transparency for overlay
negative = False        # negative conture levels
neg_color = "magenta"   # negative conture level color

first_x = False         # visibility of the first x-tick label
first_y = False         # visibility of the first x-tick label
legend = True           # visibility of the legend in overlay plots

# GRIDPLOT PARAMETERS
title_y = 0.88          # distance of the title to the axis in the grid_plot
grid_x = 0.05           # distance of the x label to the axis in the grid_plot
grid_y = 0.06           # distance of the y label to the axis in the grid_plot


# %% SPECTRA TYPE PARAMETERS

if spectra_type == "15N":
    xlim, ylim = (11.2, 6.2), (134, 104)
    x_label, y_label = "$^1$H [ppm]", "$^1$$^5$N\n[ppm]"
    y_label_grid = "$^1$$^5$N [ppm]"
    xticks, xminorticks = 1, 0.5
    yticks, yminorticks = 5, 1
    line_width = 0.3

elif spectra_type == "13C":
    xlim, ylim = (1.51, -0.58), (28, 18.2)
    x_label, y_label = "$^1$H [ppm]", "$^1$$^3$C\n[ppm]"
    y_label_grid = "$^1$$^3$C [ppm]"
    xticks, xminorticks = 0.5, 0.1
    yticks, yminorticks = 5, 1
    line_width = 0.6

elif spectra_type == "CON":
    xlim, ylim = (178, 168), (140, 108)
    x_label, y_label = "$^1$$^3$C [ppm]", "$^1$$^5$N\n[ppm]"
    y_label_grid = "$^1$$^5$N [ppm]"
    xticks, xminorticks = 2, 0.5
    yticks, yminorticks = 5, 1
    line_width = 0.3

elif spectra_type == "ZOOM":
    xlim, ylim = (8.8, 8), (118, 114)
    x_label, y_label = "$^1$H [ppm]", "$^1$$^5$N\n[ppm]"
    y_label_grid = "$^1$$^5$N [ppm]"
    xticks, xminorticks = 0.2, 0.05
    yticks, yminorticks = 2, 0.5
    line_width = 0.6

elif spectra_type == "SIZE":
    xlim, ylim = (8.9, 6.6), (130, 110)
    x_label, y_label = "$^1$H [ppm]", "$^1$$^5$N\n[ppm]"
    y_label_grid = "$^1$$^5$N [ppm]"
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

tick_width = 1
plt.rcParams["xtick.major.width"] = tick_width
plt.rcParams["ytick.major.width"] = tick_width
plt.rcParams["xtick.minor.width"] = tick_width
plt.rcParams["ytick.minor.width"] = tick_width

labels = True
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
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["mathtext.fontset"] = "dejavusans"
else:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
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
    "labels": labels,
    "alpha": alpha,
    "negative": negative,
    "neg_color": neg_color,
    "save_svg": save_svg,
    "save_png": save_png,
    "first_x": first_x,
    "first_y": first_y,
    "legend": legend
}


# %% FILES TO OVERLAY
# this gives you the number of the spectrum for the overlay plot

print("\n----- Spectrum Menu -----")
for i, name in enumerate(file_names):
    print(f" [{i}] -> {name}")
print("-------------------------\n")

# %% PLOT ALL SPECTRA
# -- you can define a "folder" for all plots --
# -- if no folder is defined everything will be saved in the "results" folder --
if True:
    pf.plot_everything(p, folder="results/single")

# %% PLOT INDIVIDUAL OVERLAYS
# -- you need to define a "name" for the overlay plot --
# -- you can plot all combinations, the last number will be on top --

if True:
    pf.overlay(p, [0, 1], name="overlay 1", folder="results/over")
    pf.overlay(p, [2, 3], name="overlay 2", folder="results/over")

# %% PLOT ALL SPECTRA IN A GRID
# -- in "grid_plot" you need to define a row and a col for the grid --
if True:
    pf.grid_plot(p, row=2, col=2, folder="results/grid")

    # %% -- in "grid_plot_over" you need to define a spectra to overlay all others --

    pf.grid_plot_over(p, over=0, row=2, col=2, folder="results/grid", reverse=False)

    # %% -- in "grid_plot_over_xp" you need to define overlay groups for each subplot --

    pf.grid_plot_over_xp(
        p,
        row=2,
        col=2,
        folder="results/grid",
        overlay_groups=[
            [0, 1],
            [2, 3],
            [1, 2, 3],
            [0, 2, 3],
        ],
    )
