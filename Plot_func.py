import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.ticker as mticker
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from adjustText import adjust_text
import os
import contextlib

# Global dictionary to store pre-calculated contour vertices
_CONTOUR_CACHE = {}


# --- HELPER FUNCTIONS ---


def parse_plot_data(data, base_dir, default_sub_dir="pdata/1"):
    """
    Parses the input data list into individual Series for plotting.
    Handles optional columns like csv_file dynamically.
    Automatically appends the default pdata folder if none is specified in the path.
    """
    df = pd.DataFrame(data)

    # Rename the first 4 columns, and grab the 5th if it exists
    if len(df.columns) >= 5:
        df.columns = ["folder", "name", "contour", "color", "csv_file"] + list(df.columns[5:])
    else:
        df.columns = ["folder", "name", "contour", "color"]
        df["csv_file"] = None  # Create an empty column if no CSVs are provided

    # If the user already wrote "pdata" in the folder name, use it as-is.
    # Otherwise, append the default_sub_dir (pdata/1).
    df["full_path"] = df["folder"].apply(
        lambda x: str(base_dir / x) if "pdata" in str(x) else str(base_dir / x / default_sub_dir)
    )

    # Return the extracted columns as a tuple
    return (df["full_path"], df["name"], df["contour"], df["color"], df["csv_file"])


def save_and_clear(fig, folder, name, p):
    """Handles multi-format saving, showing, and memory cleanup."""

    if "spectra_type" in p:
        name = f"{name}_{p['spectra_type']}"

    save_base = Path(folder) / name

    # Save PNG and SVG
    if p.get("save_png", False):
        fig.savefig(
            save_base.with_suffix(".png"),
            transparent=False,
            bbox_inches="tight",
            dpi=p.get("dpi", 300),
        )
    if p.get("save_svg", False):
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
        if len(xticklabels) >= 2:
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


def hide_texts_outside_axes(ax, texts, arrows=None):
    """
    Hide text labels whose final adjusted bounding box lies outside the axes.
    Also hides corresponding arrows if provided.
    """
    fig = ax.figure
    fig.canvas.draw()

    renderer = fig.canvas.get_renderer()
    ax_bbox = ax.get_window_extent(renderer)

    for i, text in enumerate(texts):
        text_bbox = text.get_window_extent(renderer)

        inside = (
            text_bbox.x0 >= ax_bbox.x0
            and text_bbox.x1 <= ax_bbox.x1
            and text_bbox.y0 >= ax_bbox.y0
            and text_bbox.y1 <= ax_bbox.y1
        )

        if not inside:
            text.set_visible(False)

            if arrows is not None and i < len(arrows):
                arrows[i].set_visible(False)


def add_labels_from_csv(ax, csv_file, p):
    """
    Reads positions from a CSV, creates text objects on the axis,
    and returns the objects and coordinates to be adjusted later.
    """
    if not p.get("peak_labels", True):
        return [], [], []

    if not csv_file or pd.isna(csv_file):
        return [], [], []

    csv_path = Path(p.get("csv_dir")) / csv_file

    if not csv_path.exists():
        print(f"  -> Warning: CSV file not found at {csv_path.resolve()}")
        return [], [], []

    try:
        df_labels = pd.read_csv(csv_path)

        required_cols = ["Residue", "position_1", "position_2"]
        if not all(col in df_labels.columns for col in required_cols):
            print(f"  -> Warning: CSV file {csv_file} is missing required columns.")
            return [], [], []

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        xmin, xmax = sorted(xlim)
        ymin, ymax = sorted(ylim)

        df_labels = df_labels[
            (df_labels["position_1"] >= xmin)
            & (df_labels["position_1"] <= xmax)
            & (df_labels["position_2"] >= ymin)
            & (df_labels["position_2"] <= ymax)
        ]

        if df_labels.empty:
            return [], [], []

        texts = []
        x = df_labels["position_1"].to_numpy()
        y = df_labels["position_2"].to_numpy()

        for _, row in df_labels.iterrows():
            t = ax.text(
                row["position_1"],
                row["position_2"],
                str(row["Residue"]),
                fontsize=p.get("label_fontsize", 6),
                color="black",
                ha="center",
                va="center",
                clip_on=False,
            )
            texts.append(t)

        return texts, x.tolist(), y.tolist()

    except Exception as e:
        print(f"  -> Error reading or plotting {csv_file}: {e}")
        return [], [], []


def adjust_all_labels(ax, texts, x, y, p):
    """Runs adjustText on a pooled collection of text objects from all overlaid spectra."""
    if not texts:
        return

    ax.figure.canvas.draw()
    np.random.seed(p.get("peak_seed", 42))

    with open(os.devnull, "w") as fnull:
        with contextlib.redirect_stdout(fnull):
            adjusted_texts, arrows = adjust_text(
                texts,
                x=x,
                y=y,
                ax=ax,
                prevent_crossings=True,
                ensure_inside_axes=True,
                expand_axes=False,
                expand=p.get("expand", (1.75, 1.5)),
                force_text=(0.2, 0.1),
                force_static=(0.2, 0.1),
                force_pull=(0.01, 0.01),
                iter_lim=p.get("iter_lim", 100),
                min_arrow_len=0,
                arrowprops=dict(
                    arrowstyle="-",
                    color="black",
                    lw=0.5,
                    alpha=0.5,
                    shrinkA=0.5,
                    shrinkB=0,
                ),
            )

    # Hide labels that ended up outside the axis after adjust_text
    hide_texts_outside_axes(ax, adjusted_texts, arrows)


# --- MAIN PLOTTING FUNCTIONS ---


def plot_everything(p):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    files = p["files"]
    folder = p.get("out_single")
    result_folder(folder)

    for i, f in enumerate(files):
        if dic_all[i] is None:
            continue

        fig = plt.figure(dpi=p["dpi"], figsize=(p["xsize"], p["ysize"]))
        ax = fig.add_subplot()

        apply_formatting(ax, p, title=p["file_names"][i], add_labels=True)
        draw_contours(ax, dic_all[i], data_all[i], p, p["cont"][i], p["colors"][i])

        # 2-Step Labeling
        t, x, y = add_labels_from_csv(ax, p["csv_files"][i], p)
        adjust_all_labels(ax, t, x, y, p)

        save_and_clear(fig, folder, f"{p['file_names'][i]}", p)


def overlay(p, get, name):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    folder = p.get("out_overlay")
    result_folder(folder)

    fig = plt.figure(dpi=p["dpi"], figsize=(p["xsize"], p["ysize"]))
    ax = fig.add_subplot()

    apply_formatting(ax, p, title=name, add_labels=True)

    h_all, names = [], []
    all_texts, all_x, all_y = [], [], []

    for idx in get:
        if dic_all[idx] is None:
            continue
        names.append(p["file_names"][idx])
        clp = draw_contours(ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx])

        # Collect labels from this spectrum
        t, x, y = add_labels_from_csv(ax, p["csv_files"][idx], p)
        all_texts.extend(t)
        all_x.extend(x)
        all_y.extend(y)

        h, _ = clp.legend_elements()
        if h:
            h_all.append(h[0])

    # Adjust all pooled labels together to prevent overlaps
    adjust_all_labels(ax, all_texts, all_x, all_y, p)

    if h_all and p.get("legend", True):
        leg = ax.legend(h_all, names, loc="upper left", framealpha=0.8, handlelength=1.5)
        for line in leg.get_lines():
            line.set_linewidth(2.0)

    save_and_clear(fig, folder, f"overlay_{name}", p)


def grid_plot(p, row=2, col=2):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    folder = p.get("out_grid")
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

            # 2-Step Labeling
            t, x, y = add_labels_from_csv(ax, p["csv_files"][i], p)
            adjust_all_labels(ax, t, x, y, p)
        else:
            ax.axis("off")

    save_and_clear(fig, folder, "grid", p)


def grid_plot_over(p, over, row=2, col=2, reverse=False):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    folder = p.get("out_grid")
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

            all_texts, all_x, all_y = [], [], []

            if not reverse:
                # Draw Main spectrum (Bottom)
                draw_contours(
                    ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx]
                )
                # Draw Overlay spectrum (Top)
                draw_contours(
                    ax, dic_all[over], data_all[over], p, p["cont"][over], p["colors"][over]
                )
            else:
                # Draw Overlay spectrum (Bottom)
                draw_contours(
                    ax, dic_all[over], data_all[over], p, p["cont"][over], p["colors"][over]
                )
                # Draw Main spectrum (Top)
                draw_contours(
                    ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx]
                )

            # 1. Labels for the main spectrum (idx)
            t1, x1, y1 = add_labels_from_csv(ax, p["csv_files"][idx], p)
            all_texts.extend(t1)
            all_x.extend(x1)
            all_y.extend(y1)

            # 2. Labels for the overlay spectrum (over)
            t2, x2, y2 = add_labels_from_csv(ax, p["csv_files"][over], p)
            all_texts.extend(t2)
            all_x.extend(x2)
            all_y.extend(y2)

            # Adjust labels for this subplot together
            adjust_all_labels(ax, all_texts, all_x, all_y, p)
        else:
            ax.axis("off")

    save_and_clear(fig, folder, "grid_over", p)


def grid_plot_over_xp(p, overlay_groups, row=2, col=2):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    folder = p.get("out_grid")
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
            all_texts, all_x, all_y = [], [], []

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

                # Collect labels
                t, x, y = add_labels_from_csv(ax, p["csv_files"][spec_idx], p)
                all_texts.extend(t)
                all_x.extend(x)
                all_y.extend(y)

                h, _ = clp.legend_elements()
                if h:
                    h_all.append(h[0])

            # Adjust all pooled labels together
            adjust_all_labels(ax, all_texts, all_x, all_y, p)

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
