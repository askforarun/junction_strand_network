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
  Optional helper script for one-command export/plot workflows.

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

## Installation

```bash
pip install -r requirements.txt
```

Dependencies:

- `networkx`
- `numpy`
- `matplotlib`

## Quickstart

Run from repo root (`/users/ass2009/sharedscratch/junction_strand_network`).

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

![Full junction-strand network](junction_strand_network.png)

Zoomed local junction-strand motif:

![Zoomed junction-strand motif](junction_strand_network_zoom.png)

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
| `pva_atom_u` | The PVA endpoint atom index bonded to the `source` GLU junction. |
| `pva_atom_v` | The PVA endpoint atom index bonded to the `target` GLU junction. |
| `glu_atom_u` | The reactive GLU atom index bonded to `pva_atom_u`. |
| `glu_atom_v` | The reactive GLU atom index bonded to `pva_atom_v`. |
| `is_self_loop` | `0` if the PVA strand connects two different GLU junctions; `1` if both ends connect to the same GLU junction. For the validated networks in this repo, this should be `0` for all edges. |

Legacy note:

- Older CSV files may also include `pva_chain_id`.
- In the older exporter, `pva_chain_id` was the identifier used for that PVA strand, which in practice was the lower atom index of the two reactive PVA endpoint atoms.

## Topology Checks Reported

`junction_strand_network_report.txt` includes:

- Number of junction nodes and strand edges
- Connected components and connectivity status
- Junction functionality / projected degree distribution
- 4-regular and 4-graph-connectivity checks
- Self-loop and parallel-edge diagnostics

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
