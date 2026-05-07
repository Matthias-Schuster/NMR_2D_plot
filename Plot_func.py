import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.ticker as mticker
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

# Global dictionary to store pre-calculated contour vertices
_CONTOUR_CACHE = {}


# --- HELPER FUNCTIONS ---


def save_and_clear(fig, folder, name, p):
    """Handles multi-format saving, showing, and memory cleanup."""
    save_base = Path(folder) / name
    # Save PNG and SVG
    if p["save_png"]:
        fig.savefig(
            save_base.with_suffix(".png"),
            transparent=False,
            bbox_inches="tight",
            dpi=p.get("dpi", 300),
        )
    if p["save_svg"]:
        fig.savefig(save_base.with_suffix(".svg"), transparent=False, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def apply_formatting(ax, p, title=None, is_grid=False, add_labels=False):
    """Handles titles, limits, ticks, and optionally adds axis labels."""

    if p.get("titles", True) and title:
        ax.set_title(title, y=p["title_y"] if is_grid else 1.0)
    ax.set_xlim(p["xlim"])
    ax.set_ylim(p["ylim"])
    ax.xaxis.set_major_locator(mticker.MultipleLocator(p["xticks"]))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(p["xminorticks"]))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(p["yticks"]))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(p["yminorticks"]))

    if not is_grid:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

    # Hide the first x-tick label
    if p.get("first_x", False):
        xticklabels = ax.get_xticklabels()
        if len(xticklabels) >= 2 and not is_grid:
            xticklabels[-2].set_visible(False)

    # Hide the top y-tick label
    if p.get("first_y", False):
        yticklabels = ax.get_yticklabels()
        if len(yticklabels) > 0:
            yticklabels[1].set_visible(False)

    if add_labels and p["labels"]:
        ax.set_xlabel(p["x_label"], loc="left", labelpad=p["labelpad_x"], backgroundcolor="w")
        ax.set_ylabel(
            p["y_label"],
            loc="top",
            labelpad=p["labelpad_y"],
            rotation="horizontal",
            ma="center",
            backgroundcolor="w",
        )


def draw_contours(ax, dic, data, p, contour_start, color):
    """Handles ppm scale calculation and both positive/negative contour plotting."""
    if isinstance(color, list):
        color = tuple(color)

    # 1. Create cache key
    cache_key = (id(data), contour_start, p["factor"], p["lines"], p["negative"], color)

    x_ppm, y_ppm = read_hsqc_bruker_pdata(dic, data)
    extent = (x_ppm[0], x_ppm[-1], y_ppm[0], y_ppm[-1])

    # 2. IF IN CACHE: Draw the lines directly from memory
    if cache_key in _CONTOUR_CACHE:
        cached = _CONTOUR_CACHE[cache_key]

        # Draw Positive Contours
        lc_pos = LineCollection(
            cached["pos_paths"], colors=color, linewidths=p["line_width"], alpha=p["alpha"]
        )
        ax.add_collection(lc_pos)

        # Draw Negative Contours
        if p["negative"] and cached["neg_paths"]:
            lc_neg = LineCollection(
                cached["neg_paths"],
                colors=p["neg_color"],
                linewidths=p["line_width"],
                alpha=p["alpha"],
                linestyles="dashed",
            )
            ax.add_collection(lc_neg)

        # Attach the legend dummy method directly to the object
        lc_pos.legend_elements = lambda: ([Line2D([0], [0], color=color, lw=p["line_width"])], [""])

        return lc_pos

    # 3. IF NOT IN CACHE: Perform standard calculation
    levels = contour_start * p["factor"] ** np.arange(p["lines"])

    # Plot Positive Contours
    clp = ax.contour(
        data, levels, colors=color, linewidths=p["line_width"], extent=extent, alpha=p["alpha"]
    )

    # Extract positive paths directly (Matplotlib 3.8+ syntax)
    pos_paths = []
    for path in clp.get_paths():
        pos_paths.extend(path.to_polygons())

    negative = p["negative"]
    neg_paths = []

    # Plot Negative Contours
    if negative:
        neg_levels = np.sort(-levels)
        cln = ax.contour(
            data,
            neg_levels,
            colors=p["neg_color"],
            linewidths=p["line_width"],
            extent=extent,
            alpha=p["alpha"],
        )

        # Extract negative paths directly (Matplotlib 3.8+ syntax)
        for path in cln.get_paths():
            neg_paths.extend(path.to_polygons())

    # 4. Save to global cache
    _CONTOUR_CACHE[cache_key] = {"pos_paths": pos_paths, "neg_paths": neg_paths}

    return clp


# --- MAIN PLOTTING FUNCTIONS ---


def plot_everything(p, folder="results"):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    files = p["files"]
    result_folder(folder)

    for i, f in enumerate(files):
        if dic_all[i] is None:
            continue

        fig = plt.figure(dpi=p["dpi"], figsize=(p["xsize"], p["ysize"]))
        ax = fig.add_subplot()

        apply_formatting(ax, p, title=p["file_names"][i], add_labels=True)

        draw_contours(ax, dic_all[i], data_all[i], p, p["cont"][i], p["colors"][i])

        save_and_clear(fig, folder, f"{p['file_names'][i]}", p)


def overlay(p, get, name, folder="results"):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    result_folder(folder)

    fig = plt.figure(dpi=p["dpi"], figsize=(p["xsize"], p["ysize"]))
    ax = fig.add_subplot()

    h_all, names = [], []
    for idx in get:
        if dic_all[idx] is None:
            continue
        names.append(p["file_names"][idx])
        clp = draw_contours(ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx])
        h, _ = clp.legend_elements()
        if h:
            h_all.append(h[0])

    apply_formatting(ax, p, title=name, add_labels=True)

    if h_all and p.get("legend", True):
        leg = ax.legend(h_all, names, loc="upper left", framealpha=0.8, handlelength=1.5)
        for line in leg.get_lines():
            line.set_linewidth(2.0)

    save_and_clear(fig, folder, f"overlay_{name}", p)


def grid_plot(p, row=2, col=2, folder="results"):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    result_folder(folder)

    fig, axes = plt.subplots(
        row,
        col,
        dpi=p["dpi"],
        figsize=(p["xsize"] * col, p["ysize"] * row),
        sharey=True,
        sharex=True,
        gridspec_kw={"wspace": 0, "hspace": 0},
    )

    if p["labels"]:
        fig.supxlabel(p["x_label"], y=p["grid_x"])
        fig.supylabel(p["y_label_grid"], x=p["grid_y"])

    axes_flat = np.atleast_1d(axes).flatten()
    for i, ax in enumerate(axes_flat):
        if i < len(dic_all) and dic_all[i]:
            apply_formatting(ax, p, title=p["file_names"][i], is_grid=True)
            draw_contours(ax, dic_all[i], data_all[i], p, p["cont"][i], p["colors"][i])
        else:
            ax.axis("off")

    save_and_clear(fig, folder, "grid", p)


def grid_plot_over(p, over, row=2, col=2, folder="results", reverse=False):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    result_folder(folder)

    # Prepare data (excluding the 'over' index for the grid base)
    indices = [i for i in range(len(dic_all)) if i != over]

    fig, axes = plt.subplots(
        row,
        col,
        dpi=p["dpi"],
        figsize=(p["xsize"] * col, p["ysize"] * row),
        sharey=True,
        sharex=True,
        gridspec_kw={"wspace": 0, "hspace": 0},
    )

    if p["labels"]:
        fig.supxlabel(p["x_label"], y=p["grid_x"])
        fig.supylabel(p["y_label_grid"], x=p["grid_y"])

    axes_flat = np.atleast_1d(axes).flatten()
    for i, ax in enumerate(axes_flat):
        if i < len(indices):
            idx = indices[i]
            apply_formatting(ax, p, title=p["file_names"][idx], is_grid=True)
            if dic_all[idx] is None or dic_all[over] is None:
                continue
            if not reverse:
                # Draw Main spectrum (Bottom)
                draw_contours(ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx])
                # Draw Overlay spectrum (Top)
                draw_contours(
                    ax, dic_all[over], data_all[over], p, p["cont"][over], p["colors"][over]
                )
            if reverse:
                # Draw Overlay spectrum (Bottom)
                draw_contours(
                    ax, dic_all[over], data_all[over], p, p["cont"][over], p["colors"][over]
                )
                # Draw Main spectrum (Top)
                draw_contours(ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx])
        else:
            ax.axis("off")

    save_and_clear(fig, folder, "grid_over", p)


def grid_plot_over_xp(p, overlay_groups, row=2, col=2, folder="results"):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    result_folder(folder)

    fig, axes = plt.subplots(
        row,
        col,
        dpi=p["dpi"],
        figsize=(p["xsize"] * col, p["ysize"] * row),
        sharey=True,
        sharex=True,
        gridspec_kw={"wspace": 0, "hspace": 0},
    )

    if p["labels"]:
        fig.supxlabel(p["x_label"], y=p["grid_x"])
        fig.supylabel(p["y_label_grid"], x=p["grid_y"])

    axes_flat = np.atleast_1d(axes).flatten()
    for i, ax in enumerate(axes_flat):
        if i < len(overlay_groups):
            apply_formatting(ax, p, is_grid=True)
            h_all, subplot_names = [], []

            for spec_idx in overlay_groups[i]:
                if dic_all[spec_idx] is None:
                    continue
                subplot_names.append(p["file_names"][spec_idx])
                clp = draw_contours(
                    ax,
                    dic_all[spec_idx],
                    data_all[spec_idx],
                    p,
                    p["cont"][spec_idx],
                    p["colors"][spec_idx],
                )
                h, _ = clp.legend_elements()
                if h:
                    h_all.append(h[0])

            if h_all and p.get("legend", True):
                leg = ax.legend(
                    h_all, subplot_names, loc="upper left", framealpha=0.8, handlelength=1.5
                )
                for line in leg.get_lines():
                    line.set_linewidth(2.0)
        else:
            ax.axis("off")

    save_and_clear(fig, folder, "grid_over_xp", p)


# --- UTILITIES ---


def read_hsqc_bruker_pdata(dic, data):
    x_ppm = _ppm_axis(
        data.shape[1], dic["procs"]["OFFSET"], dic["procs"]["SW_p"], dic["procs"]["SF"]
    )
    y_ppm = _ppm_axis(
        data.shape[0], dic["proc2s"]["OFFSET"], dic["proc2s"]["SW_p"], dic["proc2s"]["SF"]
    )
    return x_ppm, y_ppm


def _ppm_axis(size, offset_ppm, sw_hz, sf_mhz):
    sw_ppm = float(sw_hz) / float(sf_mhz)
    return float(offset_ppm) - (np.arange(size, dtype=float) * (sw_ppm / size))


def read_data(files):
    data_all, dic_all = [], []
    for f in files:
        path = Path(f)
        if path.exists():
            dic, data = ng.bruker.read_pdata(str(path))
            dic_all.append(dic)
            data_all.append(data)
        else:
            print(f"Path missing: {f}")
            dic_all.append(None)
            data_all.append(None)
    return dic_all, data_all


def result_folder(name):
    Path(name).mkdir(parents=True, exist_ok=True)


def colormap(files, maps="viridis"):
    cmap = plt.get_cmap(maps)
    colors_array = cmap(np.linspace(0, 1, len(files)))
    # Convert the array to a list of tuples
    return [tuple(c) for c in colors_array]
