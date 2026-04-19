#!/usr/bin/env python3
"""
Plot a coarse Junction-Strand Network from a single exported edge CSV map.

Input files are expected to come from:
    python analysis/junction_strand_network.py --data ... --pva ... --glu ... --atom1-type ... --atom2-type ... -n ... --csv-out junction_edges.csv

Nodes are GLU junctions. Edges are PVA strands.

CLI Usage Examples:
    # 1. Basic default spring plot
    $ python analysis/plot_junction_strand_network.py --csv junction_edges.csv --out plot.png

    # 2. Kamada-Kawai layout (spatially distributes nodes relative to strand lengths)
    $ python analysis/plot_junction_strand_network.py --csv junction_edges.csv --out plot_kk.png --layout kk

    # 3. High-resolution, thick-edged publication-quality plot (40x40 canvas, 300 DPI, thick lines)
    $ python analysis/plot_junction_strand_network.py --csv junction_edges.csv --out pub_plot.png --layout kk --large

    # 4. Same as above, but with white bold node IDs visibly overlaid onto the junction rings
    $ python analysis/plot_junction_strand_network.py --csv junction_edges.csv --out pub_plot_labels.png --layout kk --large --labels
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from typing import Dict, Optional, Sequence

import networkx as nx


def read_junction_strand_network_data(csv_file: str) -> nx.MultiGraph:
    """Read exported a single network edge CSV file into a MultiGraph."""
    graph = nx.MultiGraph()

    with open(csv_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            source = int(row["source"])
            target = int(row["target"])
            graph.add_node(source, kind="glu_ring", ring_id=source)
            graph.add_node(target, kind="glu_ring", ring_id=target)
            graph.add_edge(
                source,
                target,
                edge_id=int(row["edge_id"]),
                pva_chain_id=int(row["pva_chain_id"]) if row["pva_chain_id"] else None,
                pva_atom_u=int(row["pva_atom_u"]) if row["pva_atom_u"] else None,
                pva_atom_v=int(row["pva_atom_v"]) if row["pva_atom_v"] else None,
                glu_atom_u=int(row["glu_atom_u"]) if row["glu_atom_u"] else None,
                glu_atom_v=int(row["glu_atom_v"]) if row["glu_atom_v"] else None,
                is_self_loop=bool(int(row["is_self_loop"])) if row["is_self_loop"] else False,
            )

    return graph


def write_junction_strand_network_plot(
    graph: nx.MultiGraph,
    output_file: str,
    layout: str = "spring",
    with_labels: bool = False,
    with_edge_labels: bool = False,
    large: bool = False,
) -> str:
    """Write a plot of the Junction-Strand Network from exported data.
    
    Args:
        graph: NetworkX MultiGraph to plot.
        output_file: Path to save the image.
        layout: String algorithm identifier ("spring" or "kk").
        with_labels: Whether to overlay GLU junction IDs onto the nodes.
        with_edge_labels: Whether to overlay PVA strand IDs onto the edges.
        large: Scale up to a massive 40x40 high-DPI canvas with thicker edges and clear node colors.
    """
    os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    projection = nx.Graph()
    projection.add_nodes_from(graph.nodes(data=True))
    projection.add_edges_from(graph.edges())

    if projection.number_of_nodes() == 0:
        raise ValueError("Cannot plot an empty Junction-Strand Network.")

    nodelist = sorted(projection.nodes())

    if layout == "spring":
        pos = nx.spring_layout(projection, seed=7)
    elif layout == "kk":
        pos = nx.kamada_kawai_layout(projection)
    else:
        raise ValueError(f"Unsupported layout '{layout}'. Use 'spring' or 'kk'.")

    if projection.number_of_nodes() > 100:
        if large:
            node_size = 1000
            edge_alpha = 0.95
            node_linewidth = 2.0
        else:
            node_size = 50
            edge_alpha = 0.75
            node_linewidth = 0.45
    else:
        node_size = 260 if projection.number_of_nodes() <= 80 else 120
        edge_alpha = 0.9
        node_linewidth = 0.8

    fig_size = (40, 40) if large else (8, 8)
    fig, ax = plt.subplots(figsize=fig_size)
    if large:
        edge_widths = [6.0 + 2.0 * (graph.number_of_edges(u, v) - 1) for u, v in projection.edges()]
        edge_color = "#1E90FF"  # Vibrant Dodger Blue
    else:
        edge_widths = [0.8 + 0.5 * (graph.number_of_edges(u, v) - 1) for u, v in projection.edges()]
        edge_color = "#2f6f9f"

    nx.draw_networkx_edges(
        projection,
        pos,
        width=edge_widths,
        edge_color=edge_color,
        alpha=edge_alpha,
        ax=ax,
    )

    node_color = "#FF4500" if large else "#d8e8c8"  # OrangeRed when large
    node_edge_color = "#4B0082" if large else "#1f2a1f"  # Indigo border when large

    nx.draw_networkx_nodes(
        projection,
        pos,
        nodelist=nodelist,
        node_size=node_size,
        node_color=node_color,
        edgecolors=node_edge_color,
        linewidths=node_linewidth,
        ax=ax,
    )

    if with_labels:
        labels = {
            node: str(projection.nodes[node].get("ring_id", node))
            for node in nodelist
        }
        if large:
            nx.draw_networkx_labels(
                projection,
                pos,
                labels=labels,
                font_size=13,
                font_color="white",
                font_weight="bold",
                ax=ax,
            )
        else:
            nx.draw_networkx_labels(projection, pos, labels=labels, font_size=8, ax=ax)

    if with_edge_labels:
        edge_labels = {}
        for u, v in projection.edges():
            edge_data = graph.get_edge_data(u, v, default={})
            pva_ids = sorted(
                {
                    attrs.get("pva_chain_id")
                    for attrs in edge_data.values()
                    if attrs.get("pva_chain_id") is not None
                }
            )
            if pva_ids:
                edge_labels[(u, v)] = ",".join(str(pva_id) for pva_id in pva_ids)

        nx.draw_networkx_edge_labels(
            projection,
            pos,
            edge_labels=edge_labels,
            font_size=8 if large else 5,
            font_color="#3a2a00",
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.75,
                "pad": 0.15,
            },
            ax=ax,
        )

    ax.set_title(
        "Junction-Strand Network (junction nodes, strand labels)",
        fontsize=28 if large else 12,
        fontweight="bold",
    )
    ax.set_axis_off()
    fig.tight_layout()
    dpi_val = 300 if large else 220
    fig.savefig(output_file, dpi=dpi_val, bbox_inches="tight")
    plt.close(fig)
    return output_file


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint for plotting exported junction-strand network data."""
    parser = argparse.ArgumentParser(
        description="Plot a Junction-Strand Network from a single exported edge CSV map."
    )
    parser.add_argument("--csv", required=True, help="Edge CSV file from junction_strand_network.py")
    parser.add_argument("--out", required=True, help="Output plot path")
    parser.add_argument(
        "--layout",
        choices=("spring", "kk"),
        default="spring",
        help="Layout to use for plotting (default: spring). 'kk' refers to Kamada-Kawai layout.",
    )
    parser.add_argument(
        "--large",
        action="store_true",
        help="Use a very large figure size (30x30) to expand overly dense graphs.",
    )
    parser.add_argument(
        "--labels",
        action="store_true",
        help="Draw node labels (GLU junction ids).",
    )
    parser.add_argument(
        "--edge-labels",
        action="store_true",
        help="Draw edge labels (PVA strand ids).",
    )
    args = parser.parse_args(argv)

    graph = read_junction_strand_network_data(args.csv)
    write_junction_strand_network_plot(
        graph,
        args.out,
        layout=args.layout,
        with_labels=args.labels,
        with_edge_labels=args.edge_labels,
        large=args.large,
    )
    print(f"Wrote plot: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
