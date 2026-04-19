#!/usr/bin/env python3
"""
Build a coarse Junction-Strand Network from a crosslinked LAMMPS data file.

Nodes are GLU rings (junctions). Edges are PVA chains (strands) inferred from
formed crosslinks whose atom types are provided explicitly. A NetworkX MultiGraph is used so multiple distinct PVA
When exported to CSV, the generated file is an edge-list containing the following information:
    - edge_id: A unique identifier for the edge (strand).
    - source: The node_id (ring_id) of the origin GLU junction.
    - target: The node_id (ring_id) of the destination GLU junction.
    - pva_chain_id: The ID of the PVA chain connecting the two junctions.
    - pva_atom_u: The specific PVA atom index at the source endpoint.
    - pva_atom_v: The specific PVA atom index at the target endpoint.
    - glu_atom_u: The specific GLU atom index matching the source PVA endpoint.
    - glu_atom_v: The specific GLU atom index matching the target PVA endpoint.
    - is_self_loop: 1 if the source and target junctions are identical, else 0.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from collections import defaultdict
from typing import Any, Dict, Optional, Sequence, Tuple

import networkx as nx
import numpy as np

from src.system_constants import PVA_ATOMS_PER_MONOMER

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

JunctionNode = Tuple[str, int]


def read_atom_types_from_lammps_data(datafile: str) -> dict:
    """Parse the LAMMPS Atoms section as ``atom id -> atom type``."""
    atom_types = {}
    with open(datafile) as f:
        for line in f:
            if "Atoms" in line:
                for atom_line in f:
                    if atom_line.strip() and (atom_line[0].isalpha() or atom_line[0] == "_"):
                        break
                    parts = atom_line.split()
                    if len(parts) >= 3:
                        atom_types[parts[0]] = parts[2]
                break
    if not atom_types:
        raise RuntimeError(f"Failed to read atom types from Atoms section in {datafile}")
    return atom_types


def read_crosslink_edges_for_types(
    datafile: str,
    atom_types: dict,
    atom1_type: int,
    atom2_type: int,
) -> list[tuple[int, int]]:
    """Parse Bonds and return bonds between the configured crosslink atom types."""
    atom1_type = str(atom1_type)
    atom2_type = str(atom2_type)
    crosslink_edges = []
    with open(datafile) as f:
        for line in f:
            if "Bonds" in line:
                for bond_line in f:
                    if "Angles" in bond_line:
                        break
                    parts = bond_line.split()
                    if len(parts) < 4:
                        continue
                    a, b = parts[2], parts[3]
                    if a in atom_types and b in atom_types:
                        if (
                            (atom_types[a] == atom1_type and atom_types[b] == atom2_type)
                            or (atom_types[a] == atom2_type and atom_types[b] == atom1_type)
                        ):
                            crosslink_edges.append((int(a), int(b)))
                break
    return crosslink_edges


def build_pva_edges_and_chain_map(pva_file: str, shift_index: int):
    """Build PVA strand endpoint edges from ``PVA.txt``.

    Every second sorted atom id is used as the first endpoint. The second
    endpoint is offset by ``shift_index``.
    """
    pva_atoms = np.loadtxt(pva_file).astype(int)
    pva_atoms.sort()
    first_endpoints = pva_atoms[::2]
    pva_edges = [(int(u), int(u) + shift_index) for u in first_endpoints]
    all_pva_atoms = {atom for edge in pva_edges for atom in edge}
    pva_atom_to_chain = {}
    for u, v in pva_edges:
        base = min(u, v)
        pva_atom_to_chain[u] = base
        pva_atom_to_chain[v] = base
    return pva_edges, all_pva_atoms, pva_atom_to_chain


def build_glu_ring_maps(glu_file: str):
    """Build GLU ring membership maps from ``GLU.txt``."""
    glu_atoms = np.loadtxt(glu_file).astype(int)
    glu_atoms.sort()
    ring_seeds = glu_atoms.tolist()[::4]

    glu_edges = []
    for seed in ring_seeds:
        j = int(seed)
        glu_edges.extend([(j, j + 2), (j + 2, j + 12), (j + 12, j + 10), (j + 10, j)])

    glu_graph = nx.Graph()
    glu_graph.add_edges_from(glu_edges)

    glu_ring_id = {}
    ring_to_atoms = {}
    for ring_id, component in enumerate(nx.connected_components(glu_graph)):
        component_set = {int(atom) for atom in component}
        ring_to_atoms[ring_id] = component_set
        for atom in component_set:
            glu_ring_id[atom] = ring_id

    return glu_edges, glu_ring_id, ring_to_atoms, set(glu_ring_id.keys())


def build_junction_strand_network(
    data_file: str,
    pva_file: str,
    glu_file: str,
    n: int,
    atom1_type: int,
    atom2_type: int,
    allow_loop_edges: bool = False,
) -> nx.MultiGraph:
    """Build a Junction-Strand Network from formed crosslinks between configured atom types.

    Args:
        data_file: Crosslinked LAMMPS data file.
        pva_file: PVA index file used to define chain endpoints.
        glu_file: GLU index file used to define ring membership.
        n: Number of repeat units per PVA chain.
        atom1_type: First crosslink atom type in the LAMMPS data file.
        atom2_type: Second crosslink atom type in the LAMMPS data file.
        allow_loop_edges: If True, allow a PVA chain whose two endpoints map to
            the same GLU ring to be represented as a self-loop edge. If False,
            such a chain raises ValueError.

    Returns:
        A MultiGraph with GLU ring nodes and one PVA-chain edge per strand.

    Raises:
        ValueError: If a formed bond between the configured atom types
            cannot be mapped to the expected PVA endpoint / GLU ring sets, or if a participating chain does not have
            exactly one formed crosslink on each endpoint. Distinct GLU rings are
            also required unless ``allow_loop_edges`` is True.
    """
    shift_index = (PVA_ATOMS_PER_MONOMER * n - 2) - 1
    atom_types = read_atom_types_from_lammps_data(data_file)
    crosslink_edges = read_crosslink_edges_for_types(
        data_file,
        atom_types,
        atom1_type=atom1_type,
        atom2_type=atom2_type,
    )
    pva_edges, all_pva_atoms, pva_atom_to_chain = build_pva_edges_and_chain_map(
        pva_file, shift_index
    )
    _, glu_ring_id, _, _ = build_glu_ring_maps(glu_file)

    chain_endpoints = {min(u, v): (u, v) for u, v in pva_edges}
    chain_records: Dict[int, list[Tuple[int, int, int]]] = defaultdict(list)

    for a, b in crosslink_edges:
        if a in all_pva_atoms and b in glu_ring_id:
            pva_atom, glu_atom = a, b
        elif b in all_pva_atoms and a in glu_ring_id:
            pva_atom, glu_atom = b, a
        else:
            raise ValueError(
                f"Crosslink bond ({a}, {b}) could not be mapped to a PVA endpoint "
                "and a GLU ring from the provided PVA/GLU index files."
            )

        chain_id = pva_atom_to_chain[pva_atom]
        chain_records[chain_id].append((pva_atom, glu_atom, glu_ring_id[glu_atom]))

    graph = nx.MultiGraph()

    for chain_id, records in sorted(chain_records.items()):
        if chain_id not in chain_endpoints:
            raise ValueError(f"PVA chain id {chain_id} is missing from the endpoint map.")

        endpoint_u, endpoint_v = chain_endpoints[chain_id]
        by_endpoint: Dict[int, list[Tuple[int, int]]] = defaultdict(list)
        for pva_atom, glu_atom, ring_id in records:
            by_endpoint[pva_atom].append((glu_atom, ring_id))

        actual_endpoints = set(by_endpoint)
        expected_endpoints = {endpoint_u, endpoint_v}
        if actual_endpoints != expected_endpoints:
            raise ValueError(
                f"PVA chain {chain_id} must have exactly one formed "
                f"{atom1_type}-{atom2_type} bond on "
                f"each endpoint {sorted(expected_endpoints)}, found endpoints "
                f"{sorted(actual_endpoints)}."
            )

        if len(by_endpoint[endpoint_u]) != 1 or len(by_endpoint[endpoint_v]) != 1:
            raise ValueError(
                f"PVA chain {chain_id} must have exactly one formed "
                f"{atom1_type}-{atom2_type} bond per "
                "endpoint."
            )

        glu_atom_u, ring_u = by_endpoint[endpoint_u][0]
        glu_atom_v, ring_v = by_endpoint[endpoint_v][0]
        if ring_u == ring_v and not allow_loop_edges:
            raise ValueError(
                f"PVA chain {chain_id} connects twice to GLU ring {ring_u}; "
                "distinct junctions are required."
            )

        node_u: JunctionNode = ("glu_ring", ring_u)
        node_v: JunctionNode = ("glu_ring", ring_v)
        graph.add_node(node_u, kind="glu_ring", ring_id=ring_u)
        graph.add_node(node_v, kind="glu_ring", ring_id=ring_v)
        graph.add_edge(
            node_u,
            node_v,
            pva_chain_id=chain_id,
            pva_atom_u=endpoint_u,
            pva_atom_v=endpoint_v,
            glu_ring_u=ring_u,
            glu_ring_v=ring_v,
            glu_atom_u=glu_atom_u,
            glu_atom_v=glu_atom_v,
            inferred_from=((endpoint_u, glu_atom_u), (endpoint_v, glu_atom_v)),
        )

    return graph


def summarize_junction_strand_network(graph: nx.MultiGraph) -> Dict[str, Any]:
    """Summarize connectivity and path metrics for a Junction-Strand Network."""
    projection = nx.Graph()
    projection.add_nodes_from(graph.nodes(data=True))
    projection.add_edges_from(graph.edges())

    components = sorted(nx.connected_components(projection), key=len, reverse=True)
    n_components = len(components)
    largest_size = len(components[0]) if components else 0
    is_connected = n_components == 1 and projection.number_of_nodes() > 0

    n_junction_nodes = graph.number_of_nodes()
    n_strand_edges = graph.number_of_edges()
    n_unique_junction_connections = projection.number_of_edges()
    self_loop_count = sum(1 for u, v, _k in graph.edges(keys=True) if u == v)

    degree_by_node = dict(graph.degree())
    functionality_distribution: Dict[int, int] = {}
    for degree in degree_by_node.values():
        functionality_distribution[degree] = functionality_distribution.get(degree, 0) + 1
    average_functionality = (
        sum(degree_by_node.values()) / n_junction_nodes if n_junction_nodes else None
    )
    max_functionality = max(degree_by_node.values()) if degree_by_node else None

    projected_degree_by_node = dict(projection.degree())
    projected_degree_distribution: Dict[int, int] = {}
    for degree in projected_degree_by_node.values():
        projected_degree_distribution[degree] = (
            projected_degree_distribution.get(degree, 0) + 1
        )
    average_projected_degree = (
        sum(projected_degree_by_node.values()) / n_junction_nodes
        if n_junction_nodes
        else None
    )
    max_projected_degree = (
        max(projected_degree_by_node.values()) if projected_degree_by_node else None
    )
    is_projected_four_regular = (
        n_junction_nodes > 0
        and len(projected_degree_distribution) == 1
        and 4 in projected_degree_distribution
    )
    has_4_graph_connectivity = (
        is_connected
        and is_projected_four_regular
        and self_loop_count == 0
        and (n_strand_edges == n_unique_junction_connections)
    )

    avg_path: Optional[float]
    diameter: Optional[int]
    path_metric_scope: Optional[str]

    if not components:
        avg_path = None
        diameter = None
        path_metric_scope = None
    else:
        metric_graph = (
            projection if is_connected else projection.subgraph(components[0]).copy()
        )
        avg_path = nx.average_shortest_path_length(metric_graph)
        diameter = nx.diameter(metric_graph)
        path_metric_scope = "full_graph" if is_connected else "largest_component"

    return {
        "n_junction_nodes": n_junction_nodes,
        "n_strand_edges": n_strand_edges,
        "n_unique_junction_connections": n_unique_junction_connections,
        "parallel_edge_count": n_strand_edges - n_unique_junction_connections,
        "self_loop_strand_count": self_loop_count,
        "average_junction_functionality": average_functionality,
        "max_junction_functionality": max_functionality,
        "functionality_distribution": dict(sorted(functionality_distribution.items())),
        "average_projected_junction_degree": average_projected_degree,
        "max_projected_junction_degree": max_projected_degree,
        "projected_degree_distribution": dict(
            sorted(projected_degree_distribution.items())
        ),
        "is_projected_four_regular": is_projected_four_regular,
        "has_4_graph_connectivity": has_4_graph_connectivity,
        "connected_components": n_components,
        "largest_connected_component_size": largest_size,
        "is_connected": is_connected,
        "path_metric_scope": path_metric_scope,
        "average_shortest_path_length": avg_path,
        "diameter": diameter,
    }


def summarize_loop_defects(graph: nx.MultiGraph) -> Dict[str, Any]:
    """Summarize self-loop strand defects in a Junction-Strand Network.

    A self-loop strand is a PVA chain whose two endpoints are crosslinked to
    the same GLU junction. The input graph should be built with
    ``allow_loop_edges=True`` if loop defects are possible.
    """
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    node_comp = {node: cid for cid, comp in enumerate(components) for node in comp}

    component_summaries = []
    for cid, comp in enumerate(components):
        subgraph = graph.subgraph(comp)
        strand_count = subgraph.number_of_edges()
        loop_count = sum(1 for u, v, _key in subgraph.edges(keys=True) if u == v)
        component_summaries.append(
            {
                "cid": cid,
                "junctions": len(comp),
                "strands": strand_count,
                "loops": loop_count,
            }
        )

    loop_details = []
    loops_in_disconnected_components = 0
    for u, v, _key, attrs in graph.edges(keys=True, data=True):
        if u != v:
            continue
        cid = node_comp.get(u)
        ring_id = attrs["glu_ring_u"]
        loop_details.append(
            {
                "pva_atom_u": attrs["pva_atom_u"],
                "pva_atom_v": attrs["pva_atom_v"],
                "shared_rings": [ring_id],
                "component_id": cid if cid is not None else "?",
                "pva_chain_id": attrs["pva_chain_id"],
            }
        )
        if cid not in (None, 0):
            loops_in_disconnected_components += 1

    largest_component = component_summaries[0] if component_summaries else None
    largest_component_ties = []
    if largest_component:
        largest_key = (largest_component["junctions"], largest_component["strands"])
        largest_component_ties = [
            summary["cid"]
            for summary in component_summaries
            if (summary["junctions"], summary["strands"]) == largest_key
        ]

    return {
        "total_self_loop_strands": len(loop_details),
        "connected_components": len(components),
        "largest_component": largest_component,
        "largest_component_ties": largest_component_ties,
        "loops_in_disconnected_components": loops_in_disconnected_components,
        "component_summaries": component_summaries,
        "loop_details": loop_details,
    }


def format_loop_defect_summary(summary: Dict[str, Any]) -> str:
    """Return loop-defect report text for topology reports and CLI output."""
    lines = [
        "",
        "================ LOOP SUMMARY ================",
        f"Total self-loop strands detected : {summary['total_self_loop_strands']}",
        f"Junction-strand connected components: {summary['connected_components']}",
    ]

    largest = summary["largest_component"]
    if largest:
        lines.append(
            "Largest junction-strand component "
            f"(cid={largest['cid']}, junctions={largest['junctions']}, "
            f"strands={largest['strands']}, "
            f"loops={largest['loops']})"
        )
        tied_cids = summary["largest_component_ties"]
        if len(tied_cids) > 1:
            lines.append(f"Largest component tie: cids={tied_cids}")
    else:
        lines.append("Largest junction-strand component (none)")

    lines.append(
        f"Loops in disconnected components: {summary['loops_in_disconnected_components']}"
    )

    component_summaries = summary["component_summaries"]
    if component_summaries:
        lines.append("Component summaries:")
        for component_summary in component_summaries:
            lines.append(
                f"  cid={component_summary['cid']} "
                f"junctions={component_summary['junctions']} "
                f"strands={component_summary['strands']} "
                f"loops={component_summary['loops']}"
            )
    lines.extend(["================================================", ""])

    if not summary["loop_details"]:
        lines.append("No loops detected.")
    else:
        lines.append("chain_u  chain_v  shared_rings  comp_id  pva_chain_id")
        for detail in summary["loop_details"]:
            lines.append(
                f"{detail['pva_atom_u']}  {detail['pva_atom_v']}  "
                f"{detail['shared_rings']}  {detail['component_id']}  "
                f"{detail['pva_chain_id']}"
            )

    return "\n".join(lines)


def export_junction_strand_network_data(
    graph: nx.MultiGraph,
    csv_file: str,
) -> str:
    """Write an edge list table for downstream plotting or analysis.

    Columns:
    - edge_id: A unique identifier for the edge (strand).
    - source: The node_id (ring_id) of the origin GLU junction.
    - target: The node_id (ring_id) of the destination GLU junction.
    - pva_chain_id: The ID of the PVA chain connecting the two junctions.
    - pva_atom_u: The specific PVA atom index at the source endpoint.
    - pva_atom_v: The specific PVA atom index at the target endpoint.
    - glu_atom_u: The specific GLU atom index matching the source PVA endpoint.
    - glu_atom_v: The specific GLU atom index matching the target PVA endpoint.
    - is_self_loop: 1 if the source and target junctions are identical, else 0.
    """
    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "edge_id",
                "source",
                "target",
                "pva_chain_id",
                "pva_atom_u",
                "pva_atom_v",
                "glu_atom_u",
                "glu_atom_v",
                "is_self_loop",
            ],
        )
        writer.writeheader()
        for edge_id, (u, v, _key, attrs) in enumerate(
            sorted(
                graph.edges(keys=True, data=True),
                key=lambda item: (item[0][1], item[1][1], item[3].get("pva_chain_id", -1)),
            ),
            start=1,
        ):
            writer.writerow(
                {
                    "edge_id": edge_id,
                    "source": u[1],
                    "target": v[1],
                    "pva_chain_id": attrs.get("pva_chain_id", ""),
                    "pva_atom_u": attrs.get("pva_atom_u", ""),
                    "pva_atom_v": attrs.get("pva_atom_v", ""),
                    "glu_atom_u": attrs.get("glu_atom_u", ""),
                    "glu_atom_v": attrs.get("glu_atom_v", ""),
                    "is_self_loop": int(u == v),
                }
            )

    return csv_file


def _format_summary(summary: Dict[str, Any]) -> str:
    """Return terminal-friendly text for a Junction-Strand Network summary."""
    avg_path = summary["average_shortest_path_length"]
    avg_path_text = "None" if avg_path is None else f"{avg_path:.6f}"
    diameter = summary["diameter"]
    diameter_text = "None" if diameter is None else str(diameter)
    scope = summary["path_metric_scope"] or "n/a"
    avg_func = summary["average_junction_functionality"]
    avg_func_text = "None" if avg_func is None else f"{avg_func:.6f}"
    max_func = summary["max_junction_functionality"]
    max_func_text = "None" if max_func is None else str(max_func)
    functionality_distribution = ", ".join(
        f"{functionality}:{count}"
        for functionality, count in summary["functionality_distribution"].items()
    ) or "None"
    avg_proj_degree = summary["average_projected_junction_degree"]
    avg_proj_degree_text = "None" if avg_proj_degree is None else f"{avg_proj_degree:.6f}"
    max_proj_degree = summary["max_projected_junction_degree"]
    max_proj_degree_text = "None" if max_proj_degree is None else str(max_proj_degree)
    projected_degree_distribution = ", ".join(
        f"{degree}:{count}"
        for degree, count in summary["projected_degree_distribution"].items()
    ) or "None"

    return "\n".join(
        [
            "Junction-Strand Network Summary",
            "================================",
            f"GLU junction nodes: {summary['n_junction_nodes']}",
            f"PVA strand edges: {summary['n_strand_edges']}",
            f"Unique junction connections: {summary['n_unique_junction_connections']}",
            f"Parallel strand edges: {summary['parallel_edge_count']}",
            f"Self-loop strands: {summary['self_loop_strand_count']}",
            f"Average junction functionality: {avg_func_text}",
            f"Max junction functionality: {max_func_text}",
            f"Functionality distribution: {functionality_distribution}",
            f"Average projected junction degree: {avg_proj_degree_text}",
            f"Max projected junction degree: {max_proj_degree_text}",
            f"Projected degree distribution: {projected_degree_distribution}",
            f"Connected components: {summary['connected_components']}",
            f"Largest connected component size: {summary['largest_connected_component_size']}",
            f"Connected: {'YES' if summary['is_connected'] else 'NO'}",
            f"Projected graph 4-regular: {'YES' if summary['is_projected_four_regular'] else 'NO'}",
            f"4-graph connectivity: {'YES' if summary['has_4_graph_connectivity'] else 'NO'}",
            f"Path metric scope: {scope}",
            f"Average shortest path length: {avg_path_text}",
            f"Diameter: {diameter_text}",
        ]
    )


def write_junction_strand_network_report(
    summary: Dict[str, Any],
    report_file: str,
) -> str:
    """Write a text report of Junction-Strand Network properties."""
    report_path = Path(report_file)
    report_path.write_text(_format_summary(summary) + "\n")
    return str(report_path)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint for Junction-Strand Network analysis."""
    parser = argparse.ArgumentParser(
        description=(
            "Build and summarize the Junction-Strand Network "
            "(GLU junction nodes, PVA strand edges)."
        )
    )
    parser.add_argument("--data", required=True, help="Crosslinked LAMMPS data file")
    parser.add_argument("--pva", required=True, help="PVA index file (e.g. PVA.txt)")
    parser.add_argument("--glu", required=True, help="GLU index file (e.g. GLU.txt)")
    parser.add_argument("--atom1-type", type=int, required=True, help="First crosslink atom type in the data file")
    parser.add_argument("--atom2-type", type=int, required=True, help="Second crosslink atom type in the data file")
    parser.add_argument("-n", type=int, required=True, help="PVA repeat units per chain")
    parser.add_argument("--csv-out", help="Optional CSV output path for network edge data.")
    parser.add_argument("--report-out", help="Optional text report output path.")
    args = parser.parse_args(argv)

    graph = build_junction_strand_network(
        args.data,
        args.pva,
        args.glu,
        args.n,
        atom1_type=args.atom1_type,
        atom2_type=args.atom2_type,
    )
    summary = summarize_junction_strand_network(graph)
    print(_format_summary(summary))
    if args.csv_out:
        export_junction_strand_network_data(graph, args.csv_out)
        print(f"Wrote network data: {args.csv_out}")
    if args.report_out:
        write_junction_strand_network_report(summary, args.report_out)
        print(f"Wrote network report: {args.report_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
