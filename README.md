# Junction-Strand Network

Small utilities for building and plotting a coarse junction-strand graph from a
crosslinked polymer network.

The graph representation is designed for networks such as an `A2 + B4` system:

- `B4` molecule: junction node
- `A2` molecule or polymer segment: strand edge
- node degree: number of strands attached to a junction
- edge label: strand or chain ID connecting two junctions

For the PVA-GLU hydrogel workflow, this maps to:

- `GLU` ring: junction node
- `PVA` chain: strand edge
- ideal GLU functionality: four connected PVA strands

## Files

- `junction_strand_network.py`: builds the NetworkX graph, writes a CSV edge
  table, and reports connectivity/functionality checks.
- `plot_junction_strand_network.py`: plots the CSV/graph with junction IDs on
  nodes and strand IDs on edges.
- `src/system_constants.py`: small PVA/GLU constants needed by the current
  PVA-specific parser.

## Install

Create or activate a Python environment with:

```bash
pip install -r requirements.txt
```

The required packages are `networkx`, `numpy`, and `matplotlib`.

## Build The CSV

Example for the PVA-GLU hydrogel files:

```bash
python junction_strand_network.py \
  --data data.crosslinked_updated \
  --pva PVA.txt \
  --glu GLU.txt \
  --atom1-type 9 \
  --atom2-type 10 \
  -n 17 \
  --csv-out junction_edges.csv \
  --report-out junction_strand_network_report.txt
```

The output CSV columns are:

- `edge_id`: row-level strand edge ID
- `source`: source junction ID
- `target`: target junction ID
- `pva_chain_id`: strand/PVA chain ID
- `pva_atom_u`: first PVA endpoint atom
- `pva_atom_v`: second PVA endpoint atom
- `glu_atom_u`: GLU atom bonded to the first PVA endpoint
- `glu_atom_v`: GLU atom bonded to the second PVA endpoint
- `is_self_loop`: `1` if the strand connects back to the same junction

## Plot The Network

```bash
python plot_junction_strand_network.py \
  --csv junction_edges.csv \
  --out junction_strand_network.png \
  --layout kk \
  --large \
  --labels \
  --edge-labels
```

The plot is a topology verification diagram, not a molecular snapshot. It uses
an abstract graph layout so that each junction's branches can be inspected and
cross-checked against `junction_edges.csv`.

## What To Check

For an ideal loop-free `A2 + B4` network:

- each `B4`/GLU junction has degree 4
- each `A2`/PVA strand connects two distinct junctions
- there are no self-loop strands
- there are no unexpected disconnected components
- parallel strand edges are visible in the CSV and graph summary

