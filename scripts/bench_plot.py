#!/usr/bin/env python3
"""Scatter plots of wall time against peak memory from a scripts/bench.sh CSV.

    python3 scripts/bench_plot.py book/src/benchmarks.csv out.png --dataset ERR3534523
    python3 scripts/bench_plot.py book/src/benchmarks.csv out.png

With --dataset, one panel per workload shows each program's seconds and
megabytes. Without it, every dataset is drawn relative to fastp, so sets of
different sizes share the axes. --rate turns the --dataset axes into megabases
read per second and megabases read per megabyte of peak memory, so better is
up and to the right. --name VERSION=TEXT renames a version in the
legend, for example --name 0.3.0+13331c6=my-branch. Needs matplotlib.
"""
import argparse
import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e7e6e2", "#fcfcfb"
REFERENCE = "#8a8983"
# First three categorical slots: the most that stay apart in a scatter.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
# Height of a two-line point label, in pixels.
LABEL_HEIGHT = 44
# Rough width of one label character, in pixels.
LABEL_CHAR = 8.5
# Ratios marked on the axes of the relative plot.
TICKS = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
WORKLOADS = ["single-end", "paired-end", "paired-correct-merge"]


def load(path, workers):
    """(dataset, workload) -> {(program, version): row} for one worker count."""
    groups = defaultdict(dict)
    for row in csv.DictReader(open(path)):
        if int(row["workers"]) != workers or not row.get("peak_mb"):
            continue
        row["mean"], row["peak_mb"] = float(row["mean"]), float(row["peak_mb"])
        row["megabases"] = float(row.get("input_bases") or 0) / 1e6
        groups[row["dataset"], row["workload"]][row["program"], row["version"]] = row
    return groups


def release(version):
    return tuple(int(n) for n in version.split("+")[0].split("."))


def versions(groups):
    """fasterp versions, releases in order and builds with a commit last."""
    found = {v for g in groups.values() for (p, v) in g if p == "fasterp"}
    return sorted(found, key=lambda v: ("+" in v, release(v), v))


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(length=0, pad=6, colors=MUTED)


def dot(ax, x, y, color, marker):
    ax.scatter([x], [y], s=110, color=color, marker=marker, edgecolors=SURFACE, linewidths=2, zorder=4)


def legend(fig, fastp, names, colors):
    handles = [Line2D([], [], marker="D", ls="", color=REFERENCE, markersize=8, label=fastp)]
    handles += [Line2D([], [], marker="o", ls="", color=c, markersize=9, label=n) for n, c in zip(names, colors)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.055, 0.86), ncol=len(handles),
               frameon=False, fontsize=12, handletextpad=0.3, columnspacing=1.8, labelcolor=INK)


def seconds_and_megabytes(row):
    return row["mean"], row["peak_mb"], f"{row['mean']:.2f} s · {row['peak_mb']:.0f} MB"


def rates(row):
    """Megabases per second, and megabases per megabyte of peak memory."""
    speed, economy = row["megabases"] / row["mean"], row["megabases"] / row["peak_mb"]
    return speed, economy, f"{speed:,.0f} Mb/s · {economy:.2f} Mb/MB"


def one_dataset(groups, dataset, order, names, fig, axes, measure, xlabel, ylabel):
    for ax, workload in zip(axes, WORKLOADS):
        rows = groups.get((dataset, workload), {})
        style(ax)
        ax.set_title(workload, loc="left", fontsize=13, fontweight="bold", color=INK, pad=10)
        points = [(k, r) for k, r in rows.items() if k[0] == "fastp"]
        points += [(("fasterp", v), rows["fasterp", v]) for v in order if ("fasterp", v) in rows]
        points = [(k, measure(r)) for k, r in points]
        if not points:
            continue
        path = [(x, y) for (p, _), (x, y, _) in points if p == "fasterp"]
        ax.plot(*zip(*path), color=GRID, linewidth=2, zorder=2)
        right = max(x for _, (x, _, _) in points) * 1.25
        top = max(y for _, (_, y, _) in points) * 1.2
        ax.set_xlim(0, right)
        ax.set_ylim(0, top)
        labels = []
        for (program, version), (x, y, values) in points:
            fastp = program == "fastp"
            color = REFERENCE if fastp else SERIES[order.index(version)]
            dot(ax, x, y, color, "D" if fastp else "o")
            name = f"fastp {version}" if fastp else names[order.index(version)]
            px, py = ax.transData.transform((x, y))
            # the label goes left of the point when the right would leave the panel
            labels.append([py, px, x > right * 0.6, f"{name}\n{values}"])
        # Labels on the same side that would overlap form a stack: one column
        # beside the outermost point, centred on the points, with leader lines.
        boxes = []
        for side in (False, True):
            column = sorted((l for l in labels if l[2] == side), key=lambda l: l[0])
            stacks = []
            for label in column:
                if stacks and label[0] - stacks[-1][-1][0] < LABEL_HEIGHT:
                    stacks[-1].append(label)
                else:
                    stacks.append([label])
            for stack in stacks:
                edge = (min if side else max)(l[1] for l in stack)
                middle = sum(l[0] for l in stack) / len(stack)
                for i, (py, px, left, text) in enumerate(stack):
                    x = edge + (-12 if left else 12)
                    y = py if len(stack) == 1 else middle + (i - (len(stack) - 1) / 2) * LABEL_HEIGHT
                    width = max(len(line) for line in text.split("\n")) * LABEL_CHAR
                    span = (x - width, x) if left else (x, x + width)
                    # a label that would cross one already placed takes the nearest free row
                    rows = [y + step * LABEL_HEIGHT / 2 for n in range(8) for step in (n, -n)]
                    y = next((row for row in rows if not any(
                        span[0] < hi and lo < span[1] and abs(row - other) < LABEL_HEIGHT for lo, hi, other in boxes)), y)
                    boxes.append((span[0], span[1], y))
                    ax.annotate(text, ax.transData.inverted().transform((px, py)),
                                xytext=((x - px) * 72 / fig.dpi, (y - py) * 72 / fig.dpi), textcoords="offset points",
                                fontsize=10.5, color=INK, va="center", ha="right" if left else "left",
                                linespacing=1.35, zorder=5,
                                arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.8, "shrinkB": 6}
                                if y != py else None)
        ax.set_xlabel(xlabel, color=MUTED, fontsize=11)
    axes[0].set_ylabel(ylabel, color=MUTED, fontsize=11)


def relative(groups, order, fig, axes):
    datasets = sorted({d for d, _ in groups})
    ratios = [1.0]
    for rows in groups.values():
        base = next((r for (p, _), r in rows.items() if p == "fastp"), None)
        for (p, v), r in rows.items():
            if base and p == "fasterp" and v in order:
                ratios += [r["mean"] / base["mean"], r["peak_mb"] / base["peak_mb"]]
    limits = (min(ratios) / 1.3, max(ratios) * 1.3)
    for ax, workload in zip(axes, WORKLOADS):
        style(ax)
        ax.set_title(workload, loc="left", fontsize=13, fontweight="bold", color=INK, pad=10)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(1, color=REFERENCE, linewidth=1, zorder=1)
        ax.axhline(1, color=REFERENCE, linewidth=1, zorder=1)
        for dataset in datasets:
            rows = groups.get((dataset, workload), {})
            base = next((r for (p, _), r in rows.items() if p == "fastp"), None)
            if not base:
                continue
            path = [(rows["fasterp", v]["mean"] / base["mean"], rows["fasterp", v]["peak_mb"] / base["peak_mb"], v)
                    for v in order if ("fasterp", v) in rows]
            ax.plot([x for x, _, _ in path], [y for _, y, _ in path], color=GRID, linewidth=1.5, zorder=2)
            for x, y, v in path:
                dot(ax, x, y, SERIES[order.index(v)], "o")
        dot(ax, 1, 1, REFERENCE, "D")
        ax.set_xticks(TICKS, [f"{t:g}x" for t in TICKS])
        ax.set_yticks(TICKS, [f"{t:g}x" for t in TICKS])
        ax.minorticks_off()
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_xlabel("Wall time relative to fastp", color=MUTED, fontsize=11)
    axes[0].set_ylabel("Peak memory relative to fastp", color=MUTED, fontsize=11)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("csv")
    parser.add_argument("out")
    parser.add_argument("--dataset")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--rate", action="store_true", help="plot megabases per second and per megabyte")
    parser.add_argument("--name", action="append", default=[], metavar="VERSION=TEXT")
    args = parser.parse_args()

    groups = load(args.csv, args.workers)
    # a version keeps its colour whichever datasets are drawn
    order = versions(groups)[-len(SERIES):]
    if args.dataset:
        groups = {k: g for k, g in groups.items() if k[0] == args.dataset}
    renamed = dict(n.split("=", 1) for n in args.name)
    names = [renamed.get(v, f"fasterp {v}") for v in order]
    drawn = [v in versions(groups) for v in order]
    fastp = sorted({f"fastp {v}" for g in groups.values() for (p, v) in g if p == "fastp"})
    machine = sorted({r["machine"] for g in groups.values() for r in g.values()})

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 11})
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.6), sharey=bool(not args.dataset), facecolor=SURFACE,
                             dpi=110)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.74, bottom=0.1, wspace=0.16 if args.dataset else 0.07)
    better = "lower left is better"
    if args.dataset and args.rate:
        one_dataset(groups, args.dataset, order, names, fig, axes, rates,
                    "Megabases read per second", "Megabases read per MB of peak memory")
        title = f"Throughput and memory economy on {args.dataset}"
        note = "Each line joins the fasterp builds from oldest to newest"
        better = "upper right is better"
    elif args.dataset:
        one_dataset(groups, args.dataset, order, names, fig, axes, seconds_and_megabytes,
                    "Wall time (seconds)", "Peak memory (MB)")
        title = f"Wall time and peak memory on {args.dataset}"
        note = "Each line joins the fasterp builds from oldest to newest"
    else:
        relative(groups, order, fig, axes)
        title = "Wall time and peak memory relative to fastp"
        note = f"One point per dataset ({len({d for d, _ in groups})} datasets), joined from the oldest fasterp build to the newest"
    fig.text(0.055, 0.965, title, fontsize=20, fontweight="bold", color=INK, va="top")
    fig.text(0.055, 0.9, f"gzip in, gzip out · -w {args.workers} · {' · '.join(machine)} · {note} · "
             f"{better}", fontsize=11.5, color=MUTED, va="top")
    legend(fig, " / ".join(fastp), [n for n, d in zip(names, drawn) if d],
           [c for c, d in zip(SERIES, drawn) if d])
    fig.savefig(args.out, facecolor=SURFACE)


if __name__ == "__main__":
    main()
