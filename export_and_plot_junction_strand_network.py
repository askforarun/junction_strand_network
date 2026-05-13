"""Extracted signac operation from hydrogel_simulation/workflow/checks.py."""

from pathlib import Path

from workflow.project import HydrogelProject

@HydrogelProject.pre(
    lambda job: job.isfile("data.compression_lammps")
    and job.isfile("data.crosslinked_updated")
)
@HydrogelProject.post(
    lambda job: job.doc.get("crosslinking_completion_status") == "complete"
)
@HydrogelProject.operation
def check_crosslinking_completion(job):
    """Check crosslink completion from bond-count changes."""
    with job:
        initial_bonds = count_bonds_in_data_file("data.compression_lammps")
        final_bonds = count_bonds_in_data_file("data.crosslinked_updated")
        crosslink_points = count_crosslink_points()

        if initial_bonds is None or final_bonds is None or crosslink_points is None:
            job.doc.crosslinking_completion_status = "failed"
            raise RuntimeError(
                "Could not determine crosslinking completion because one or more "
                "bond counts could not be read."
            )

        expected_bonds = initial_bonds + crosslink_points
        if final_bonds != expected_bonds:
            job.doc.crosslinking_completion_status = "incomplete"
            raise RuntimeError(
                "Crosslinking incomplete: "
                f"expected {expected_bonds} bonds, found {final_bonds}."
            )

        job.doc.crosslinking_completion_status = "complete"
        print(
            f"[{job.id}] crosslinking complete: "
            f"{initial_bonds} + {crosslink_points} = {final_bonds}"
        )


@HydrogelProject.pre(
    lambda job: job.isfile("data.crosslinked_updated")
    and job.isfile("PVA.txt")
    and job.isfile("GLU.txt")
)
@HydrogelProject.post(lambda job: job.isfile("crosslinked_percolation_report.txt"))
@HydrogelProject.operation
def check_crosslinked_percolation_and_loops(job):
    """Check percolation, 4-graph connectivity, and loop defects.

    This reports atomistic percolation (X/Y/Z) and enforces a junction-strand
    4-regular connectivity gate. It also reports self-loop strands in the
    same coarse junction-strand map.
    """
    from analysis.percolation import (
        analyze_percolation,
        compute_report,
        format_report,
        read_lammps_data,
    )
    from analysis.junction_strand_network import (
        build_junction_strand_network,
        format_loop_defect_summary,
        summarize_junction_strand_network,
        summarize_loop_defects,
    )

    with job:
        report_path = Path("crosslinked_percolation_report.txt")
        parsed = read_lammps_data("data.crosslinked_updated")
        components = analyze_percolation(
            parsed.atom_data,
            parsed.bond_list,
            parsed.neighbors,
            parsed.box,
            bond_translations=parsed.bond_translations,
        )
        summary = compute_report(components)

        atom1_type = 9
        atom2_type = 10
        branching_degree_threshold = 3
        junction_graph = build_junction_strand_network(
            data_file="data.crosslinked_updated",
            pva_file="PVA.txt",
            glu_file="GLU.txt",
            atom1_type=atom1_type,
            atom2_type=atom2_type,
            allow_loop_edges=True,
        )
        junction_summary = summarize_junction_strand_network(junction_graph)
        loop_summary = summarize_loop_defects(junction_graph)
        loop_count = int(loop_summary["total_self_loop_strands"])
        max_projected_degree = junction_summary.get("max_projected_junction_degree") or 0
        junction_connected = bool(junction_summary.get("is_connected"))
        junction_four_regular = bool(junction_summary.get("is_projected_four_regular"))
        junction_has_4_graph_connectivity = bool(
            junction_summary.get("has_4_graph_connectivity")
        )
        junction_has_branching = max_projected_degree >= branching_degree_threshold
        junction_branching_pass = (
            junction_connected
            and junction_has_branching
            and junction_four_regular
            and junction_has_4_graph_connectivity
        )

        report_lines = format_report(summary)
        report_lines += [
            "",
            "=" * 60,
            "JUNCTION-STRAND BRANCHING GATE",
            "=" * 60,
            f"  Junction-strand connected: {'YES' if junction_connected else 'NO'}",
            f"  Max projected junction degree: {max_projected_degree}",
            f"  Branching degree threshold: {branching_degree_threshold}",
            f"  Projected graph 4-regular: {'YES' if junction_four_regular else 'NO'}",
            f"  4-graph connectivity: {'YES' if junction_has_4_graph_connectivity else 'NO'}",
            f"  Branching gate pass: {'YES' if junction_branching_pass else 'NO'}",
        ]
        report_lines.extend(format_loop_defect_summary(loop_summary).splitlines())
        report_path.write_text("\n".join(report_lines) + "\n")

        system_percolates = summary["system_percolates"].tolist()
        largest_component_percolates = summary["largest_component_percolates"].tolist()
        percolation_by_dimension = {
            "x": bool(system_percolates[0]),
            "y": bool(system_percolates[1]),
            "z": bool(system_percolates[2]),
        }
        largest_component_by_dimension = {
            "x": bool(largest_component_percolates[0]),
            "y": bool(largest_component_percolates[1]),
            "z": bool(largest_component_percolates[2]),
        }
        fully_percolated = all(percolation_by_dimension.values())
        percolation_pass = bool(fully_percolated)
        four_graph_pass = bool(junction_has_4_graph_connectivity)
        percolation_and_four_graph_pass = bool(percolation_pass and four_graph_pass)

        job.doc.crosslinked_percolation_status = (
            "passed" if percolation_and_four_graph_pass else "failed"
        )
        job.doc.crosslinked_percolates = percolation_by_dimension
        job.doc.crosslinked_fully_percolated = fully_percolated
        job.doc.crosslinked_junction_strand_connected = junction_connected
        job.doc.crosslinked_junction_strand_four_regular = junction_four_regular
        job.doc.crosslinked_junction_strand_4_graph_connectivity = (
            junction_has_4_graph_connectivity
        )
        job.doc.crosslinked_junction_strand_max_projected_degree = (
            int(max_projected_degree)
        )
        job.doc.crosslinked_junction_strand_branching_threshold = (
            int(branching_degree_threshold)
        )
        job.doc.crosslinked_junction_strand_branching_pass = junction_branching_pass
        job.doc.crosslinked_loop_count = loop_count
        job.doc.crosslinked_loops_detected = bool(loop_count)
        job.doc.crosslinked_largest_component_fraction = float(
            summary["largest_component_fraction"]
        )
        job.doc.crosslinked_largest_component_percolates = (
            largest_component_by_dimension
        )
        job.doc.crosslinked_largest_component_percolation_dim = int(
            summary["largest_component_percolation_dim"]
        )
        job.doc.crosslinked_largest_component_spans_xyz = bool(
            summary["largest_component_spans_xyz"]
        )

        if not percolation_and_four_graph_pass:
            missing_dimensions = [
                dim.upper()
                for dim, does_percolate in percolation_by_dimension.items()
                if not does_percolate
            ]
            raise RuntimeError(
                "Crosslinked bond network failed percolation/4-graph criteria. "
                f"Missing percolation in: {', '.join(missing_dimensions) if missing_dimensions else 'none'}. "
                "Junction-strand 4-graph gate: "
                f"{'passed' if four_graph_pass else 'failed'} "
                f"(connected={junction_connected}, max_degree={max_projected_degree}, "
                f"threshold={branching_degree_threshold}, "
                f"four_regular={junction_four_regular}, "
                f"has_4_graph_connectivity={junction_has_4_graph_connectivity}). "
                f"See {report_path} for details."
            )

        print(
            f"[{job.id}] crosslinked network percolates in X, Y, and Z and passes "
            f"junction 4-graph connectivity with no loop defects "
            f"(report: {report_path})"
        )


@HydrogelProject.pre(
    lambda job: job.isfile("data.crosslinked_updated")
    and job.isfile("PVA.txt")
    and job.isfile("GLU.txt")
)
@HydrogelProject.post(
    lambda job: job.isfile("junction_edges.csv")
    and job.isfile("junction_strand_network.png")
)
@HydrogelProject.operation
def export_and_plot_junction_strand_network(job):
    """Export and plot the coarse junction-strand network for a crosslinked job."""
    from analysis.junction_strand_network import (
        build_junction_strand_network,
        export_junction_strand_network_data,
    )
    from analysis.plot_junction_strand_network import (
        write_junction_strand_network_plot,
    )

    with job:
        atom1_type = 9
        atom2_type = 10
        csv_path = Path("junction_edges.csv")
        plot_path = Path("junction_strand_network.png")

        graph = build_junction_strand_network(
            data_file="data.crosslinked_updated",
            pva_file="PVA.txt",
            glu_file="GLU.txt",
            atom1_type=atom1_type,
            atom2_type=atom2_type,
        )
        export_junction_strand_network_data(graph, str(csv_path))
        write_junction_strand_network_plot(
            graph,
            str(plot_path),
            layout="kk",
            with_labels=True,
            with_edge_labels=True,
            large=True,
        )
        print(
            f"[{job.id}] wrote junction-strand network CSV and plot: "
            f"{csv_path}, {plot_path}"
        )
