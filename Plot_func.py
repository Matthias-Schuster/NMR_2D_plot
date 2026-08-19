import warnings
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.ticker as mticker
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.transforms import Bbox
from adjustText import adjust_text
import os
import contextlib
import copy

# Plot_func V1.3
# Includes 3D Strip Overlays, Dynamic Negative Peak Coloring and Dynamic Spectra Styles


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
        fig.savefig(save_base.with_suffix(".svg"), transparent=True, bbox_inches="tight")

    if plt.get_backend().lower() != "agg":
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
        ax.set_xlabel(p["x_label"], loc="left", labelpad=p["labelpad_x"], backgroundcolor="None")
        ax.set_ylabel(
            p["y_label"],
            loc="top",
            labelpad=p["labelpad_y"],
            rotation="horizontal",
            ma="center",
            backgroundcolor="None",
        )


def draw_contours(ax, dic, data, p, contour_start, color):
    """Handles ppm scale calculation and both positive/negative contour plotting."""

    # 1. Determine if 'color' contains a pair for [positive, negative]
    is_two_colors = (
        isinstance(color, (list, tuple))
        and len(color) == 2
        and isinstance(color[0], (str, list, tuple))
    )

    if is_two_colors:
        pos_color = color[0]
        neg_color = color[1]
        negative = True
    else:
        pos_color = color
        neg_color = None
        negative = False

    # Convert lists to tuples to make them hashable for the cache key
    if isinstance(pos_color, (list, np.ndarray)):
        pos_color = tuple(pos_color)
    if isinstance(neg_color, (list, np.ndarray)):
        neg_color = tuple(neg_color)

    # 2. Create cache key
    cache_key = (id(data), contour_start, p["factor"], p["lines"], negative, pos_color, neg_color)

    x_ppm, y_ppm = read_hsqc_bruker_pdata(dic, data)
    extent = (x_ppm[0], x_ppm[-1], y_ppm[0], y_ppm[-1])

    # 3. IF IN CACHE: Draw the lines directly from memory
    if cache_key in _CONTOUR_CACHE:
        cached = _CONTOUR_CACHE[cache_key]

        # Draw Positive Contours
        lc_pos = LineCollection(
            cached["pos_paths"], colors=pos_color, linewidths=p["line_width"], alpha=p["alpha"]
        )
        ax.add_collection(lc_pos)

        # Draw Negative Contours
        if negative and cached["neg_paths"]:
            lc_neg = LineCollection(
                cached["neg_paths"],
                colors=neg_color,
                linewidths=p["line_width"],
                alpha=p["alpha"],
                linestyles="dashed",
            )
            ax.add_collection(lc_neg)

        # Attach the legend dummy method directly to the object
        lc_pos.legend_elements = lambda: (
            [Line2D([0], [0], color=pos_color, lw=p["line_width"])],
            [""],
        )

        return lc_pos

    # 4. IF NOT IN CACHE: Perform standard calculation
    levels = contour_start * p["factor"] ** np.arange(p["lines"])

    # Plot Positive Contours
    clp = ax.contour(
        data, levels, colors=pos_color, linewidths=p["line_width"], extent=extent, alpha=p["alpha"]
    )

    # Extract positive paths directly
    pos_paths = []
    for path in clp.get_paths():
        pos_paths.extend(path.to_polygons())

    neg_paths = []

    # Plot Negative Contours
    if negative:
        neg_levels = np.sort(-levels)
        cln = ax.contour(
            data,
            neg_levels,
            colors=neg_color,
            linewidths=p["line_width"],
            extent=extent,
            alpha=p["alpha"],
        )

        # Extract negative paths directly
        for path in cln.get_paths():
            neg_paths.extend(path.to_polygons())

    # 5. Save to global cache
    _CONTOUR_CACHE[cache_key] = {"pos_paths": pos_paths, "neg_paths": neg_paths}

    return clp


# --- LABEL FUNCTIONS ---


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

        for row in df_labels.itertuples():
            t = ax.text(
                row.position_1,
                row.position_2,
                str(row.Residue),
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


# --- LABEL ADJUSTMENT HELPER FUNCTIONS ---


def hide_texts_outside_axes(ax, texts, arrows=None):
    """Hide text labels and arrows whose bounding box lies outside the axes."""
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


def _ellipse_bboxes_from_data(ax, x, y, width, height, expand=1.0):
    """Compute display-space bounding boxes for ellipses directly from data coords."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    half_w, half_h = width / 2.0, height / 2.0
    corners_1 = np.column_stack([x - half_w, y - half_h])
    corners_2 = np.column_stack([x + half_w, y + half_h])

    disp_1 = ax.transData.transform(corners_1)
    disp_2 = ax.transData.transform(corners_2)

    x0 = np.minimum(disp_1[:, 0], disp_2[:, 0])
    x1 = np.maximum(disp_1[:, 0], disp_2[:, 0])
    y0 = np.minimum(disp_1[:, 1], disp_2[:, 1])
    y1 = np.maximum(disp_1[:, 1], disp_2[:, 1])

    if expand != 1.0:
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        hw, hh = 0.5 * (x1 - x0) * expand, 0.5 * (y1 - y0) * expand
        x0, x1, y0, y1 = cx - hw, cx + hw, cy - hh, cy + hh

    return [Bbox.from_extents(a, b, c, d) for a, b, c, d in zip(x0, y0, x1, y1)]


def _bbox_list_to_arrays(bboxes):
    """Convert list of Bbox objects to numpy arrays for fast vectorized checks."""
    n = len(bboxes)
    x0, y0, x1, y1 = np.empty(n), np.empty(n), np.empty(n), np.empty(n)
    for i, bbox in enumerate(bboxes):
        x0[i], y0[i], x1[i], y1[i] = bbox.x0, bbox.y0, bbox.x1, bbox.y1
    return x0, y0, x1, y1


def push_texts_out_of_obstacle_bboxes_fast(
    fig, ax, texts, obstacle_bboxes, padding_px=5, max_iter=80
):
    """Vectorized hard constraint step. Pushes labels out of peak ellipses."""
    if not texts or not obstacle_bboxes:
        return

    obs_x0, obs_y0, obs_x1, obs_y1 = _bbox_list_to_arrays(obstacle_bboxes)

    for _ in range(max_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        moved_any = False

        for text in texts:
            if not text.get_visible():
                continue

            tb = text.get_window_extent(renderer)
            tx0, tx1 = min(tb.x0, tb.x1), max(tb.x0, tb.x1)
            ty0, ty1 = min(tb.y0, tb.y1), max(tb.y0, tb.y1)

            # Vectorized overlap check against all obstacles simultaneously
            overlaps = (tx0 < obs_x1) & (tx1 > obs_x0) & (ty0 < obs_y1) & (ty1 > obs_y0)

            if not np.any(overlaps):
                continue

            # Calculate pixel moves that would separate the text from each overlapping obstacle side
            ox0, oy0, ox1, oy1 = (
                obs_x0[overlaps],
                obs_y0[overlaps],
                obs_x1[overlaps],
                obs_y1[overlaps],
            )

            move_left = ox0 - tx1 - padding_px
            move_right = ox1 - tx0 + padding_px
            move_down = oy0 - ty1 - padding_px
            move_up = oy1 - ty0 + padding_px

            candidate_moves = []
            for dx in move_left:
                candidate_moves.append((dx, 0.0))
            for dx in move_right:
                candidate_moves.append((dx, 0.0))
            for dy in move_down:
                candidate_moves.append((0.0, dy))
            for dy in move_up:
                candidate_moves.append((0.0, dy))

            # Pick the shortest escape route
            dx_px, dy_px = min(candidate_moves, key=lambda move: abs(move[0]) + abs(move[1]))

            old_x, old_y = text.get_position()
            old_display = ax.transData.transform((old_x, old_y))
            new_display = old_display + np.array([dx_px, dy_px])
            new_data = ax.transData.inverted().transform(new_display)

            text.set_position(new_data)
            moved_any = True

            if not moved_any:
                break


def draw_final_connectors(ax, x, y, texts, min_distance_px=2, lw=0.5, alpha=0.5, color="black"):
    """Draw connector lines once, after final label positions are known."""
    arrows = []
    for px, py, text in zip(x, y, texts):
        if not text.get_visible():
            arrows.append(None)
            continue

        tx, ty = text.get_position()

        # Calculate pixel distance to decide if we need a line at all
        anchor_disp = ax.transData.transform((px, py))
        text_disp = ax.transData.transform((tx, ty))
        distance_px = np.hypot(text_disp[0] - anchor_disp[0], text_disp[1] - anchor_disp[1])

        if distance_px <= min_distance_px:
            arrows.append(None)
            continue

        arrow = ax.annotate(
            "",
            xy=(px, py),
            xytext=(tx, ty),
            arrowprops=dict(
                arrowstyle="-",
                color=color,
                lw=lw,
                alpha=alpha,
                patchA=text,
                shrinkA=0.25,
                shrinkB=0,
            ),
            annotation_clip=False,
            zorder=4,
        )
        arrows.append(arrow)

    return arrows


# --- MAIN LABEL ADJUSTMENT FUNCTION ---


def adjust_all_labels(ax, texts, x, y, p):
    """Adjust peak labels while avoiding peak ellipses using iterative cooling."""
    if not texts:
        return [], []

    fig = ax.figure
    fig.canvas.draw()

    n = min(len(texts), len(x), len(y))
    if n == 0:
        return [], []

    texts = texts[:n]
    x = np.asarray(x[:n], dtype=float)
    y = np.asarray(y[:n], dtype=float)

    # User Settings
    ellipse_width = p.get("peak_ellipse_x", 0.06)
    ellipse_height = p.get("peak_ellipse_y", 0.6)
    show_ellipse = p.get("show_ellipse", False)
    random_seed = p.get("random_seed", 42)
    jitter = p.get("jitter", 0.3)

    # Internal parameters for placement
    bbox_expand = 1.05
    initial_offset_scale = 1.05
    expand = (1.25, 1.25)
    obstacle_padding_px = 1

    # Iteration Control ("Cooling" parameters)
    cycles = p.get("cycles", 10)  # Number of push/relax cycles
    base_max_move = 30.0  # Starting movement limit in pixels
    base_pull_force = 0.1  # Starting pull back to the peak
    iterations = p.get("iterations", 10)  # Number of max iterations

    # 1. Build obstacle boxes
    obstacle_bboxes = _ellipse_bboxes_from_data(
        ax=ax, x=x, y=y, width=ellipse_width, height=ellipse_height, expand=bbox_expand
    )

    if show_ellipse:
        for px, py in zip(x, y):
            ellipse = Ellipse(
                (px, py),
                width=ellipse_width,
                height=ellipse_height,
                fill=True,
                facecolor="black",
                edgecolor="black",
                alpha=0.25,
                zorder=1,
            )
            ax.add_patch(ellipse)

    # 2. Initial deterministic label placement (Spiral algorithm)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))

    # Fetch the starting angle from parameters and convert to radians
    user_shift_deg = p.get("starting_angle", 0.0)
    starting_phase = np.radians(user_shift_deg)

    # Apply the user's deterministic shift to the spiral
    angles = (np.arange(n, dtype=float) * golden_angle) + starting_phase

    # Calculate the base deterministic offsets
    base_dx = np.cos(angles) * ellipse_width * initial_offset_scale
    base_dy = np.sin(angles) * ellipse_height * initial_offset_scale

    # Generate reproducible random noise based on your seed
    rng = np.random.default_rng(random_seed)

    # Create random shifts scaled by the ellipse size and your chosen jitter strength
    jitter_x = rng.uniform(-jitter, jitter, n) * ellipse_width
    jitter_y = rng.uniform(-jitter, jitter, n) * ellipse_height

    # Combine the rotated spiral with the random jitter
    dx = base_dx + jitter_x
    dy = base_dy + jitter_y

    for i, text in enumerate(texts):
        text.set_position((x[i] + dx[i], y[i] + dy[i]))

    fig.canvas.draw()

    # 3. Iterative Solver Loop
    for cycle in range(cycles):

        # A. Hard constraint: Push entirely out of peaks
        push_texts_out_of_obstacle_bboxes_fast(
            fig=fig,
            ax=ax,
            texts=texts,
            obstacle_bboxes=obstacle_bboxes,
            padding_px=obstacle_padding_px,
            max_iter=iterations,
        )

        # B. Soft constraint: Resolve text overlaps
        # Decrease movement allowance and pull force as cycles progress
        current_max_move = base_max_move / (cycle + 1)
        current_pull = base_pull_force / (cycle + 1) if cycle < (cycles - 1) else 0.0

        with open(os.devnull, "w") as fnull, contextlib.redirect_stdout(
            fnull
        ), warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            adjust_text(
                texts,
                target_x=x,
                target_y=y,
                avoid_self=False,
                ax=ax,
                prevent_crossings=True,
                ensure_inside_axes=True,
                expand_axes=False,
                expand=expand,
                force_text=(0.6, 0.6),
                force_static=(0.0, 0.0),
                force_pull=(current_pull, current_pull),
                max_move=(current_max_move, current_max_move),
                iter_lim=iterations,
            )

    # 4. One final hard push to guarantee no peak is covered after the final adjust_text
    push_texts_out_of_obstacle_bboxes_fast(
        fig=fig,
        ax=ax,
        texts=texts,
        obstacle_bboxes=obstacle_bboxes,
        padding_px=obstacle_padding_px,
        max_iter=iterations,
    )

    # 5. Draw final connector lines
    arrows = draw_final_connectors(
        ax=ax,
        x=x,
        y=y,
        texts=texts,
        min_distance_px=4,
        lw=p.get("connector_width", 0.5),
        alpha=0.5,
        color="black",
    )

    hide_texts_outside_axes(ax, texts, arrows)

    return texts, arrows


# --- MAIN PLOTTING FUNCTIONS ---


def plot_everything(p):
    dic_all = p["dic_all"]
    data_all = p["data_all"]
    files = p["files"]
    folder = p.get("out_single")
    result_folder(folder)

    for i, f in enumerate(files):
        if dic_all[i] is None or data_all[i].ndim != 2:
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
        if dic_all[idx] is None or data_all[idx].ndim != 2:
            print(f"  -> Skipping [{idx}] {p['file_names'][idx]}: not a 2D spectrum.")
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


def grid_plot(p, row=2, col=2, name="grid"):
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
        if i < len(dic_all) and dic_all[i] is not None and data_all[i].ndim == 2:
            apply_formatting(ax, p, title=p["file_names"][i], is_grid=True)
            draw_contours(ax, dic_all[i], data_all[i], p, p["cont"][i], p["colors"][i])

            # 2-Step Labeling
            t, x, y = add_labels_from_csv(ax, p["csv_files"][i], p)
            adjust_all_labels(ax, t, x, y, p)
        else:
            ax.axis("off")

    save_and_clear(fig, folder, name, p)


def grid_plot_over(p, over, row=2, col=2, reverse=False, name="grid_over"):
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
            # Skip if either the base spectrum or the overlay spectrum are missing or not 2D
            if (
                dic_all[idx] is None
                or data_all[idx].ndim != 2
                or dic_all[over] is None
                or data_all[over].ndim != 2
            ):
                continue
            if dic_all[idx] is None or dic_all[over] is None:
                continue

            all_texts, all_x, all_y = [], [], []

            if not reverse:
                # Draw Main spectrum (Bottom)
                draw_contours(ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx])
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
                draw_contours(ax, dic_all[idx], data_all[idx], p, p["cont"][idx], p["colors"][idx])

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

    save_and_clear(fig, folder, name, p)


def grid_plot_over_xp(p, overlay_groups, row=2, col=2, name="grid_over_xp"):
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
                if dic_all[spec_idx] is None or data_all[spec_idx].ndim != 2:
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

    save_and_clear(fig, folder, name, p)


# --- 3D STRIP FUNCTIONS ---


def get_dim_ppm(dic_3d, data_3d, dim):
    """Calculates the PPM axis for a 3D dimension (0: F1/proc3s, 1: F2/proc2s, 2: F3/procs)."""
    if dim == 0:
        return _ppm_axis(
            data_3d.shape[0],
            dic_3d["proc3s"]["OFFSET"],
            dic_3d["proc3s"]["SW_p"],
            dic_3d["proc3s"]["SF"],
        )
    elif dim == 1:
        return _ppm_axis(
            data_3d.shape[1],
            dic_3d["proc2s"]["OFFSET"],
            dic_3d["proc2s"]["SW_p"],
            dic_3d["proc2s"]["SF"],
        )
    elif dim == 2:
        return _ppm_axis(
            data_3d.shape[2],
            dic_3d["procs"]["OFFSET"],
            dic_3d["procs"]["SW_p"],
            dic_3d["procs"]["SF"],
        )
    else:
        raise ValueError("dim must be 0, 1, or 2")


def extract_2d_strip(dic_3d, data_3d, z_slice=None, slice_index=None, slice_axis=1):
    """
    Extracts a 2D plane from a 3D Bruker dataset by PPM or point index.
    Remaps the dictionary headers to match standard 2D processing functions.

    slice_axis:
      0: Fix F1 (proc3s) -> extracts F2-F3 plane
      1: Fix F2 (proc2s) -> extracts F1-F3 plane (most common, e.g. 15N slice in HNCACB)
      2: Fix F3 (procs)  -> extracts F1-F2 plane
    """
    ppm_axis = get_dim_ppm(dic_3d, data_3d, slice_axis)

    # If PPM is provided, find the closest point index automatically
    if slice_index is None:
        if z_slice is None:
            raise ValueError("Provide either z_slice or slice_index.")
        slice_index = int(np.argmin(np.abs(ppm_axis - z_slice)))

    actual_ppm = ppm_axis[slice_index]
    dic_2d = copy.deepcopy(dic_3d)

    if slice_axis == 0:
        data_2d = data_3d[slice_index, :, :]
    elif slice_axis == 1:
        data_2d = data_3d[:, slice_index, :]
        dic_2d["proc2s"] = dic_3d["proc3s"]  # Map F1 (proc3s) to Y-axis (proc2s)
    elif slice_axis == 2:
        data_2d = data_3d[:, :, slice_index]
        dic_2d["procs"] = dic_3d["proc2s"]
        dic_2d["proc2s"] = dic_3d["proc3s"]
    else:
        raise ValueError("slice_axis must be 0, 1, or 2")

    return dic_2d, data_2d, actual_ppm, slice_index


def plot_strip(
    p,
    indices,
    z_slice=None,
    slice_index=None,
    slice_axis=1,
    name=None,
    reverse=False,
    x_slice=None,
    strip_width=0.3,
):
    """
    Plots a 2D strip from one or more 3D spectra.
    Accepts a single integer (e.g., indices=0) or a list (e.g., indices=[0, 1]).
    Can be optionally centered on an x_slice with a specific strip_width.
    """
    # Standardize input to a list
    if isinstance(indices, int):
        indices = [indices]

    # Use the guaranteed out_strip path from setup (redundant code removed!)
    folder = p["out_strip"]
    result_folder(folder)

    slices_data = []
    actual_ppm_base = None
    base_name_parts = []

    # Extract 2D slices for all requested 3D spectra
    for idx in indices:
        dic = p["dic_all"][idx]
        data = p["data_all"][idx]

        if dic is None or data.ndim != 3:
            print(f"  -> Skipping [{idx}]: missing or not 3D.")
            continue

        dic_2d, data_2d, actual_ppm, _ = extract_2d_strip(
            dic, data, z_slice, slice_index, slice_axis
        )

        # Use the ppm of the first valid spectrum for the title
        if actual_ppm_base is None:
            actual_ppm_base = actual_ppm

        slices_data.append((idx, dic_2d, data_2d))
        base_name_parts.append(p["file_names"][idx])

    if not slices_data:
        print("  -> Error: No valid 3D spectra found.")
        return

    # --- NEW NAMING LOGIC ---
    is_overlay = len(slices_data) > 1
    style = p.get("spectra_type", "")
    style_suffix = f"_{style}" if style else ""

    if name:
        base_name = name
        save_name = f"{base_name}{style_suffix}_{actual_ppm_base:.2f}ppm"
    else:
        if is_overlay:
            base_name = f"Overlay_{'_'.join(base_name_parts)}"
            save_name = f"{base_name}{style_suffix}_strip_overlay_{actual_ppm_base:.2f}ppm"
        else:
            base_name = base_name_parts[0]
            save_name = f"{base_name}{style_suffix}_strip_{actual_ppm_base:.2f}ppm"

    title = f"{base_name} ({actual_ppm_base:.2f} ppm)"
    if x_slice is not None:
        title += f"\n{x_slice:.2f} ppm"

    p_temp = p.copy()
    if x_slice is not None:
        half_width = strip_width / 2.0
        p_temp["xlim"] = (x_slice + half_width, x_slice - half_width)

    fig = plt.figure(dpi=p_temp["dpi"], figsize=(p_temp["xsize"], p_temp["ysize"]))
    ax = fig.add_subplot()

    apply_formatting(ax, p_temp, title=title, add_labels=True)

    h_all = []
    names_order = []

    # Determine plotting order (reverse list if reverse=True)
    plot_data = reversed(slices_data) if reverse else slices_data

    # Draw all contours
    for idx, dic_2d, data_2d in plot_data:
        c = draw_contours(ax, dic_2d, data_2d, p_temp, p_temp["cont"][idx], p_temp["colors"][idx])
        names_order.append(p["file_names"][idx])

        h, _ = c.legend_elements()
        if h:
            h_all.append(h[0])

    if x_slice is not None:
        ax.axvline(x=x_slice, color="black", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)

    # Apply Legend (only if it's an overlay and legends are enabled)
    if is_overlay and h_all and p_temp.get("legend", True):
        leg = ax.legend(h_all, names_order, loc="upper left", framealpha=0.8, handlelength=1.5)
        for line in leg.get_lines():
            line.set_linewidth(2.0)

    # --- NEW SAVE CALL ---
    # We pop spectra_type out of p_temp so save_and_clear doesn't append it again at the very end
    p_temp.pop("spectra_type", None)
    save_and_clear(fig, folder, save_name, p_temp)


def strip_grid_plot(
    p,
    indices,
    slices,
    slice_axis=1,
    row=1,
    col=3,
    name="strip_grid",
    strip_width=0.3,
    reverse=False,
):
    """
    Plots a grid of 2D strips from one or more 3D spectra.
    Accepts a single integer (e.g., indices=0) or a list (e.g., indices=[0, 1]).
    """
    if isinstance(indices, int):
        indices = [indices]

    folder = p["out_strip"]
    result_folder(folder)

    valid_indices = []
    for idx in indices:
        if p["dic_all"][idx] is not None and p["data_all"][idx].ndim == 3:
            valid_indices.append(idx)
        else:
            print(f"  -> Skipping [{idx}]: missing or not 3D.")

    if not valid_indices:
        print("  -> Error: No valid 3D spectra found.")
        return

    is_overlay = len(valid_indices) > 1
    if is_overlay and name == "strip_grid":
        name = "strip_grid_overlay"

    fig, axes = plt.subplots(
        row,
        col,
        dpi=p["dpi"],
        figsize=(p["xsize"] * col, p["ysize"] * row),
        sharey=True,
        sharex=False,
        gridspec_kw={"wspace": 0, "hspace": 0},
    )

    if p["labels"]:
        fig.supxlabel(p["x_label"], y=p["grid_x"])
        fig.supylabel(p["y_label_grid"], x=p["grid_y"])

    axes_flat = np.atleast_1d(axes).flatten()

    for i, ax in enumerate(axes_flat):
        if i < len(slices):
            slice_val = slices[i]

            if isinstance(slice_val, dict):
                sl_z = slice_val.get("z_slice")
                sl_idx = slice_val.get("slice_index")
                sl_title = slice_val.get("title")
                sl_x = slice_val.get("x_slice")
            else:
                sl_z = float(slice_val)
                sl_idx = None
                sl_title = None
                sl_x = None

            extracted_data = []
            actual_ppm_base = None

            for idx in valid_indices:
                dic_2d, data_2d, actual_ppm, _ = extract_2d_strip(
                    p["dic_all"][idx], p["data_all"][idx], sl_z, sl_idx, slice_axis
                )
                if actual_ppm_base is None:
                    actual_ppm_base = actual_ppm
                extracted_data.append((idx, dic_2d, data_2d))

            # Titel-Formatierung: Zeigt Name + Z-Slice + X-Slice an
            if sl_title:
                title = f"{sl_title}\n({actual_ppm_base:.2f} ppm)"
            else:
                title = f"{actual_ppm_base:.2f} ppm"

            if sl_x is not None:
                title += f"\n{sl_x:.2f} ppm"

            p_temp = p.copy()
            if sl_x is not None:
                half_width = strip_width / 2.0
                p_temp["xlim"] = (sl_x + half_width, sl_x - half_width)

            apply_formatting(ax, p_temp, title=title, is_grid=True)

            plot_data = reversed(extracted_data) if reverse else extracted_data

            h_all = []
            names_order = []

            for idx, dic_2d, data_2d in plot_data:
                c = draw_contours(
                    ax,
                    dic_2d,
                    data_2d,
                    p_temp,
                    p_temp["cont"][idx],
                    p_temp["colors"][idx],
                )
                names_order.append(p["file_names"][idx])

                h, _ = c.legend_elements()
                if h:
                    h_all.append(h[0])

            if sl_x is not None:
                ax.axvline(
                    x=sl_x, color="black", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0
                )

            if not ax.get_subplotspec().is_last_row():
                ax.tick_params(labelbottom=False)

            # Legende nur im ersten Subplot (i == 0) bei Overlays
            if i == 0 and is_overlay and h_all and p_temp.get("legend", True):
                leg = ax.legend(
                    h_all,
                    names_order,
                    loc="upper left",
                    framealpha=0.8,
                    handlelength=1.5,
                )
                for line in leg.get_lines():
                    line.set_linewidth(2.0)

        else:
            ax.axis("off")

    save_and_clear(fig, folder, name, p)


# --- UTILITIES ---


def set_style(p_dict, style_name, styles_dict=None):
    """
    Updates the plot parameters dictionary in-place with a specific spectra style.
    Looks for styles_dict passed directly, or in p_dict['SPECTRA_STYLES'].
    """
    # 1. Determine where to get the styles from
    target_styles = styles_dict or p_dict.get("SPECTRA_STYLES")

    if not target_styles:
        print("  -> Error: No SPECTRA_STYLES dictionary found in p or passed to set_style!")
        return

    # 2. Apply the chosen style
    if style_name in target_styles:
        p_dict.update(target_styles[style_name])
        p_dict["spectra_type"] = style_name  # Keeps track of the name for file exports
    else:
        print(f"  -> Warning: Style '{style_name}' not found in SPECTRA_STYLES!")


def spectrum_menu(p):
    """Prints a terminal menu of loaded spectra and their dimensionality."""
    print("\n----- Spectrum Menu -----")
    for i, name in enumerate(p["file_names"]):
        data = p["data_all"][i]
        if data is not None:
            dim_text = f"{data.ndim}D"
        else:
            dim_text = "Missing/Error"
        print(f" [{i}] -> {name} ({dim_text})")
    print("-------------------------\n")


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
