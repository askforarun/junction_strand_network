# Junction-Strand Network (A2+B4)

Python utility for reducing a crosslinked polymer network to a coarse
**junction-strand graph**, exporting the graph as a CSV edge table, and plotting
the network for topology checks.

![Junction-strand network plot](docs/junction-strand-network.png)

## Overview

A junction-strand network represents a crosslinked material as a graph:

- **junction (node)**: a multifunctional crosslinker, branch point, or network
  node
- **strand (edge)**: a polymer segment, spacer, or bifunctional molecule that
  connects two junctions
- **node degree**: the number of strands attached to a junction
- **edge label**: the molecule, chain, or strand ID connecting two junctions

This is a topology verification diagram, not a molecular snapshot. It is useful
for checking network connectivity and functionality independently of the full
atomistic geometry.

The current implementation works for `A2 + B4` networks: two-ended strands
connected to four-functional junctions. The graph abstraction can be extended to
other network models, but those cases require updates to the input mapping and
expected functionality checks.

## Key Features

- Build a NetworkX `MultiGraph` from crosslinked LAMMPS data and reactive atom
  ID files.
- Export one CSV row per strand (edge), including junction endpoints, strand IDs,
  endpoint atoms, and self-loop flags.
- Plot the abstract network with junction labels and strand labels.
- Report topology checks for connectivity, junction functionality, self-loops,
  dangling strands, and repeated connections.
- Separate the general graph model from the current PVA-glutaraldehyde input
  mapping.

## What This Tool Does

This tool converts an atomistic crosslinked network into a graph that can be
inspected directly:

1. Read a crosslinked LAMMPS data file.
2. Read two text files that define the reactive atoms used for graph mapping:
   one file, for example `PVA.txt`, lists the strand-end atoms that participate
   in crosslinking (two atoms per strand/edge), and the other, for example
   `GLU.txt`, lists the junction atoms that participate in crosslinking
   (four atoms per junction/node).
3. Identify crosslinks between the configured strand and junction atom types.
   The strand crosslinking atom type must be different from the junction
   crosslinking atom type.
4. Map each two-ended strand to the two junctions it connects.
5. Write `junction_edges.csv`.
6. Write `junction_strand_network_report.txt`.
7. Plot `junction_strand_network.png`.

### Input to Output Mapping

- `data.crosslinked`: crosslinked LAMMPS data file containing atom types and
  bonds.
- `PVA.txt`: strand-end atom IDs that participate in crosslinking; this example
  has two reactive PVA atoms per strand (edge).
- `GLU.txt`: junction atom IDs that participate in crosslinking; this example
  has four reactive glutaraldehyde atoms per junction (node).
- The strand-end atom IDs and junction atom IDs must correspond to two
  different crosslinking atom types in the LAMMPS data file.
- `junction_edges.csv`: exported edge table with one row per strand.
- `junction_strand_network.png`: abstract graph plot for visual topology
  checks.
- `junction_strand_network_report.txt`: text report summarizing network
  topology checks.

### Atom ID Ordering

The parser sorts the atom IDs in `PVA.txt` and `GLU.txt` internally, so the
file order does not matter. The sorted IDs must still follow the expected
structure:

- `PVA.txt`: two reactive strand-end atom IDs per strand (edge). After sorting,
  the parser pairs consecutive IDs as one strand:
  `(pva_1, pva_2)`, `(pva_3, pva_4)`, `(pva_5, pva_6)`, and so on.
- `GLU.txt`: four reactive junction atom IDs per junction (node). After
  sorting, the parser groups consecutive IDs as one junction:
  `(glu_1, glu_2, glu_3, glu_4)`, `(glu_5, glu_6, glu_7, glu_8)`, and so on.

For example, if sorted `PVA.txt` begins with `1, 248, 251, 498`, the first
strand uses endpoints `(1, 248)` and the second uses `(251, 498)`.
If sorted `GLU.txt` begins with `75001, 75003, 75011, 75013, 75032, 75034,
75042, 75044`, the first junction uses `(75001, 75003, 75011, 75013)` and the
second uses `(75032, 75034, 75042, 75044)`.

## Table of Contents

- Overview
- Key Features
- What This Tool Does
- Atom ID Ordering
- What You Need
- Installation
- Current Input Mapping
- Command Reference
- Tutorial
- Example Report Output
- How It Works
- Output CSV Columns
- What To Check
- Files

## What You Need

- Python 3 with `networkx`, `numpy`, and `matplotlib`.
- A crosslinked LAMMPS data file named `data.crosslinked`.
- A strand atom ID file, such as `PVA.txt`, with two reactive atoms per
  strand (edge).
- A junction atom ID file, such as `GLU.txt`, with four reactive atoms per
  junction (node).
- Two different crosslinking atom types in the LAMMPS data file: one for the
  strand atoms and one for the junction atoms.

## Installation

Create or activate a Python environment, then install the package requirements:

```bash
pip install -r requirements.txt
```

Required Python packages:

- `networkx`
- `numpy`
- `matplotlib`

## Current Input Mapping

The graph model is general. The parser currently included in this repository
maps the parent workflow's PVA-glutaraldehyde hydrogel files into a
junction-strand graph.

The example considered here is a PVA-glutaraldehyde hydrogel. The end atoms of
each PVA chain are reactive, so PVA is treated as the bifunctional `A2` strand.
The glutaraldehyde crosslinker is tetrafunctional, so it is treated as the `B4`
junction.

For that system:

- glutaraldehyde (`GLU`): junction (node)
- `PVA` chain: strand (edge)
- ideal glutaraldehyde functionality: four connected PVA strands

To use the same graph idea for another input format or network model, update
the mapping layer so it can identify junction IDs, strand IDs, and strand
endpoint crosslinks from that system's files.

## Command Reference

### Build the CSV

| Argument | Required | Description |
| --- | --- | --- |
| `--data` | yes | Crosslinked LAMMPS data file, e.g. `data.crosslinked`. |
| `--pva` | yes | PVA crosslinking atom ID file, e.g. `PVA.txt`. |
| `--glu` | yes | Glutaraldehyde crosslinking atom ID file, e.g. `GLU.txt`. |
| `--atom1-type` | yes | First crosslink atom type in the LAMMPS data file. |
| `--atom2-type` | yes | Second, different crosslink atom type in the LAMMPS data file. |
| `--csv-out` | optional | Output CSV edge table. |
| `--report-out` | optional | Output topology report. |

### Plot the Network

| Argument | Required | Description |
| --- | --- | --- |
| `--csv` | yes | Edge CSV from `junction_strand_network.py`. |
| `--out` | yes | Output PNG path. |
| `--layout` | optional | Graph layout: `spring` or `kk`. |
| `--large` | optional | Use a large high-resolution figure. |
| `--labels` | optional | Draw junction (node) labels. |
| `--edge-labels` | optional | Draw strand (edge) labels. |

## Tutorial

The tutorial assumes you are running from this repository and have
`data.crosslinked`, `PVA.txt`, and `GLU.txt` in the current directory.

### 1. Build the Junction-Strand CSV

```bash
python junction_strand_network.py \
  --data data.crosslinked \
  --pva PVA.txt \
  --glu GLU.txt \
  --atom1-type 9 \
  --atom2-type 10 \
  --csv-out junction_edges.csv \
  --report-out junction_strand_network_report.txt
```

### Example Report Output

For the example files included in this repository, the generated topology report
is:

```text
Junction-Strand Network Summary
================================
GLU junction nodes: 150
PVA strand edges: 300
Unique junction connections: 300
Parallel strand edges: 0
Self-loop strands: 0
Average junction functionality: 4.000000
Max junction functionality: 4
Functionality distribution: 4:150
Average projected junction degree: 4.000000
Max projected junction degree: 4
Projected degree distribution: 4:150
Connected components: 1
Largest connected component size: 150
Connected: YES
Projected graph 4-regular: YES
4-graph connectivity: YES
Path metric scope: full_graph
Average shortest path length: 19.127517
Diameter: 38
```

### 2. Plot the Network

```bash
python plot_junction_strand_network.py \
  --csv junction_edges.csv \
  --out junction_strand_network.png \
  --layout kk \
  --large \
  --labels \
  --edge-labels
```

The plot uses an abstract graph layout. Node labels are junction IDs and edge
labels are strand IDs, so the figure can be cross-checked directly against the
CSV.

## How It Works

1. Parse atom types from the LAMMPS `Atoms` section.
2. Parse crosslink bonds from the LAMMPS `Bonds` section using the configured
   atom types.
3. Build strand endpoint pairs from the PVA atom IDs in `PVA.txt`.
4. Build junction IDs from the glutaraldehyde atom IDs in `GLU.txt`.
5. Assign each strand endpoint crosslink to a junction.
6. Add one graph edge for every strand connecting two junctions.
7. Export the graph as a CSV edge table.
8. Summarize connectivity, functionality, loops, and parallel edges.

## Output CSV Columns

Current PVA-glutaraldehyde CSV columns:

- `edge_id`: row-level strand (edge) ID
- `source`: source junction ID
- `target`: target junction ID
- `pva_chain_id`: strand/PVA chain ID
- `pva_atom_u`: first PVA endpoint atom
- `pva_atom_v`: second PVA endpoint atom
- `glu_atom_u`: glutaraldehyde atom bonded to the first PVA endpoint
- `glu_atom_v`: glutaraldehyde atom bonded to the second PVA endpoint
- `is_self_loop`: `1` if the strand connects back to the same junction

For a fully general version, the current system-specific column names can be
renamed to generic labels such as `strand_id`, `strand_endpoint_u`,
`strand_endpoint_v`, `junction_atom_u`, and `junction_atom_v`.

## What To Check

For an ideal loop-free `A2 + B4` network:

- every `B4` junction has degree 4
- every `A2` strand connects two distinct junctions
- there are no self-loop strands
- there are no dangling strands
- there are no unexpected disconnected components
- parallel strands (edges) are visible in the CSV and graph summary

For the PVA-glutaraldehyde hydrogel example, the corresponding checks are:

- every glutaraldehyde junction should have four PVA branches
- every PVA strand should connect two distinct glutaraldehyde junctions
- the network should be connected
- self-loop and unexpected parallel strand defects should be absent unless they
  are intentionally allowed

## Files

- `junction_strand_network.py`: builds the NetworkX graph, writes a CSV edge
  table, and reports connectivity/functionality checks.
- `plot_junction_strand_network.py`: plots a CSV/graph with junction IDs on
  nodes and strand IDs on edges.
- `requirements.txt`: minimal Python dependencies.
- `data.crosslinked`: example crosslinked LAMMPS data file.
- `PVA.txt`: example PVA crosslinking atom ID file.
- `GLU.txt`: example glutaraldehyde crosslinking atom ID file.
