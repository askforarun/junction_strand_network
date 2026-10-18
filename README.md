# Junction-Strand Network (A2+B4)

Standalone Python tools to convert a crosslinked atomistic network into a
coarse junction-strand graph, export strand-level connectivity as CSV, and
plot the network for topology verification.

## What This Repo Contains

- `junction_strand_network.py`
  Builds the coarse graph from a LAMMPS data file plus reactive atom index
  files, and writes:
  - `junction_edges.csv`
  - `junction_strand_network_report.txt`
- `plot_junction_strand_network.py`
  Reads `junction_edges.csv` and writes a network plot PNG.
- `export_and_plot_junction_strand_network.py`
  Standalone command-line helper for one-command export/plot workflows.
  It does not import signac or require the hydrogel simulation repository.

## Graph Definition

- Node (`junction`): one GLU junction (tetrafunctional crosslinker unit)
- Edge (`strand`): one PVA strand connecting two junctions
- Node labels in plot: GLU junction IDs
- Edge labels in plot: sequential `edge_id` values from the CSV

This is a topology representation, not a molecular-geometry snapshot.

## Inputs

- Crosslinked LAMMPS data file (example: `data.crosslinked`)
- `PVA.txt`: reactive PVA endpoint atom indices (2 per strand)
- `GLU.txt`: reactive GLU atom indices (4 per junction)
- Two distinct reactive atom types in the LAMMPS data

Current mapping in this repo is PVA-GLU specific (`A2+B4`).
The graph is inferred from explicit reactive indices and formed crosslinks,
not from a hard-coded atom count per repeat unit. It therefore supports the
corrected fully hydrolyzed PVA strands (one hydroxyl per two backbone carbons)
and the current 31-atom GLU model junctions. GLU denotes the model junction
building block rather than free glutaraldehyde.

Atom IDs are positive, one-based LAMMPS atom IDs. PVA IDs are sorted and
paired, while GLU IDs are sorted and grouped in fours. This requires atom
numbering in which each strand's endpoints and each junction's four sites
occupy consecutive groups in those sorted lists. Arbitrarily interleaved or
renumbered molecules require regenerated index files and membership handling.
Every strand endpoint and every listed GLU site must have exactly one formed
crosslink; incomplete strands and duplicate reactive IDs are rejected.

## Current production sizes

The root-level `data.crosslinked`, `PVA.txt`, `GLU.txt`, reports, and diagram
are a small demonstration with 150 junctions and 300 strands. They are not
the current production boxes. Each production replica has:

| Strand length | PVA strands / graph edges | GLU junctions / graph nodes | Atoms |
|---|---:|---:|---:|
| 17 | 2700 | 1350 | 371,250 |
| 21 | 2100 | 1050 | 347,550 |
| 25 | 1800 | 900 | 348,300 |

Each strand length has three independently prepared replicas. Supplied
`examples/production_n*/` directories contain self-contained CSV edge lists,
reactive index files, and topology reports for one current replica per strand
length. They can be plotted without access to the simulation repository;
re-exporting from atomic data requires the corresponding LAMMPS input.

## Installation

```bash
pip install -r requirements.txt
```

Dependencies:

- `networkx`
- `numpy`
- `matplotlib`

## Quickstart

Run from this repository's root directory.

### Combined standalone command

```bash
python export_and_plot_junction_strand_network.py \
  --data data.crosslinked --pva PVA.txt --glu GLU.txt --output-dir results
```

Defaults are reactive atom types 9 and 10 and a reproducible spring layout.
For large networks, omit labels and avoid `--large` unless the high-resolution
canvas is needed. Use `--skip-plot` to export only CSV and topology metrics.

To plot a current production example:

```bash
python plot_junction_strand_network.py \
  --csv examples/production_n25/junction_edges.csv \
  --out results_n25.png --layout spring
```

Regression tests, including incomplete-input rejection:

```bash
python -m unittest discover -s tests -v
```

### 1. Build CSV + report

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

### 2. Plot network

```bash
python plot_junction_strand_network.py \
  --csv junction_edges.csv \
  --out junction_strand_network.png \
  --layout kk \
  --large \
  --labels \
  --edge-labels
```

## Network Diagrams

Full junction-strand network:

![Full junction-strand network](docs/junction_strand_network_readme.png)

Zoomed local junction-strand motif:

![Zoomed junction-strand motif](docs/junction_strand_network_zoom_readme.png)

## CSV Schema and Column Meanings

Current `junction_edges.csv` columns are:

- `edge_id`
- `source`
- `target`
- `pva_atom_u`
- `pva_atom_v`
- `glu_atom_u`
- `glu_atom_v`
- `is_self_loop`

Detailed meanings:

| Column | Meaning |
| --- | --- |
| `edge_id` | A simple sequential edge number in the exported CSV: `1, 2, 3, ...`. This is an edge/row label. |
| `source` | The GLU junction ID at one end of the PVA strand. |
| `target` | The GLU junction ID at the other end of the PVA strand. |
| `pva_atom_u` | The PVA endpoint atom index that forms the covalent PVA-GLU C-C crosslink bond at the `source` junction side. |
| `pva_atom_v` | The PVA endpoint atom index that forms the covalent PVA-GLU C-C crosslink bond at the `target` junction side. |
| `glu_atom_u` | The reactive GLU carbon atom index covalently bonded to `pva_atom_u` (the source-side PVA-GLU C-C bond pair). |
| `glu_atom_v` | The reactive GLU carbon atom index covalently bonded to `pva_atom_v` (the target-side PVA-GLU C-C bond pair). |
| `is_self_loop` | `0` if the PVA strand connects two different GLU junctions; `1` if both ends connect to the same GLU junction. For the validated networks in this repo, this should be `0` for all edges. |

Short example from the diagram:

- Consider junction nodes `62` and `64` connected by edge label `129` in the zoomed network image.
- The corresponding CSV row is:
  `129,62,64,15617,15736,38535,38585,0`
- Interpretation:
  `edge_id = 129` is the plotted strand label, `source = 62` and `target = 64` are the two GLU junction nodes, (`pva_atom_u`, `glu_atom_u`) = (`15617`, `38535`) is one covalent PVA-GLU C-C bond pair, (`pva_atom_v`, `glu_atom_v`) = (`15736`, `38585`) is the other end bond pair, and `is_self_loop = 0` confirms this strand links two distinct junctions.

## Topology Checks Reported

`junction_strand_network_report.txt` includes:

- Number of junction nodes and strand edges
- Connected components and connectivity status
- Junction functionality / projected degree distribution
- 4-regular and 4-graph-connectivity checks
- Self-loop and parallel-edge diagnostics

The legacy `4-graph connectivity` field means connected, projected
4-regular, and free of self-loops and parallel strands. It does not compute
mathematical four-vertex connectivity. Graph connectivity also does not
establish spatial percolation across periodic boundaries. Likewise, "no
self-loop defects" does not mean absence of multi-strand graph cycles,
which are expected in a connected crosslinked network.

### Current Report Interpretation (Example in This Repo)

For the current example data shipped in this repo, the report values are:

- `GLU junction nodes: 150`
- `PVA strand edges: 300`
- `Parallel strand edges: 0`
- `Self-loop strands: 0`
- `Functionality distribution: 4:150`
- `Connected components: 1`
- `Connected: YES`
- `Projected graph 4-regular: YES`
- `4-graph connectivity: YES`

Interpretation:

- The network is fully connected as one macroscopic component.
- Every GLU junction has the expected functionality of 4 (ideal `B4` behavior).
- No self-loop or parallel-strand defects are detected in this sample.
- This is consistent with the specified connected `A2+B4` junction-strand topology.

## Typical Validation Targets (Ideal Loop-Free A2+B4)

- Junction degree = 4 for all GLU junctions
- No self-loops (`is_self_loop = 0` for all edges)
- Single connected component
- No unintended parallel strands unless explicitly allowed

## Files

- `junction_strand_network.py`
- `plot_junction_strand_network.py`
- `export_and_plot_junction_strand_network.py`
- `requirements.txt`
- `data.crosslinked` (example input)
- `PVA.txt` (example reactive strand endpoints)
- `GLU.txt` (example reactive junction atoms)
