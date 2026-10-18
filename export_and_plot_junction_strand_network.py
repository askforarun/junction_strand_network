"""Standalone command-line export, topology report, and optional network plot."""
import argparse
from pathlib import Path

from junction_strand_network import (
    build_junction_strand_network,
    export_junction_strand_network_data,
    summarize_junction_strand_network,
    write_junction_strand_network_report,
)
from plot_junction_strand_network import write_junction_strand_network_plot


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True)
    parser.add_argument("--pva", required=True)
    parser.add_argument("--glu", required=True)
    parser.add_argument("--atom1-type", type=int, default=9)
    parser.add_argument("--atom2-type", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--layout", choices=("spring", "kk"), default="spring")
    parser.add_argument("--labels", action="store_true")
    parser.add_argument("--edge-labels", action="store_true")
    parser.add_argument("--large", action="store_true")
    parser.add_argument("--skip-plot", action="store_true")
    parser.add_argument("--allow-loop-edges", action="store_true")
    args = parser.parse_args(argv)
    graph = build_junction_strand_network(
        args.data, args.pva, args.glu, args.atom1_type, args.atom2_type,
        allow_loop_edges=args.allow_loop_edges,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv = args.output_dir / "junction_edges.csv"
    report = args.output_dir / "junction_strand_network_report.txt"
    export_junction_strand_network_data(graph, str(csv))
    write_junction_strand_network_report(summarize_junction_strand_network(graph), str(report))
    print(report.read_text())
    print(f"Wrote {csv}")
    if not args.skip_plot:
        plot = args.output_dir / "junction_strand_network.png"
        write_junction_strand_network_plot(
            graph, str(plot), layout=args.layout, with_labels=args.labels,
            with_edge_labels=args.edge_labels, large=args.large,
        )
        print(f"Wrote {plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
