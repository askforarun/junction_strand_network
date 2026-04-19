# Junction-Strand Network

Utilities for representing a crosslinked polymer network as a coarse
**junction-strand graph**.

The idea is general:

- **junction node**: a multifunctional crosslinker, branch point, or network
  node
- **strand edge**: a polymer segment, spacer, or bifunctional molecule that
  connects two junctions
- **node degree**: the number of strands attached to a junction
- **edge label**: the molecule, chain, or strand ID connecting two junctions

This representation is useful for checking network topology independently of
the full atomistic geometry. It is a topology verification diagram, not a
molecular snapshot.

## Example: A2 + B4 Network

For an ideal `A2 + B4` network:

- `B4` molecule: junction node
- `A2` molecule or polymer segment: strand edge
- ideal `B4` functionality: degree 4
- each `A2` strand should connect two distinct `B4` junctions

The coarse graph therefore turns:

```text
B4_i -- A2_k -- B4_j
```

into:

```text
node i ---- edge k ---- node j
```

This makes it easy to detect:

- dangling strands
- self-loop strands
- under-functionalized junctions
- over-functionalized junctions
- disconnected components
- parallel edges between the same junction pair

## Current Parser

The graph concept is general, but the current parser included here is
specialized for the PVA-GLU hydrogel files produced by the parent workflow.
For that system:

- `GLU` ring: junction node
- `PVA` chain: strand edge
- ideal GLU functionality: four connected PVA strands

The PVA-GLU parser expects:

- a crosslinked LAMMPS data file
- `PVA.txt`, defining PVA endpoint atoms
- `GLU.txt`, defining GLU junction/ring atoms
- the two crosslinking atom types used in the LAMMPS data file

To use the same graph idea for another chemistry, only the mapping layer needs
to be changed: define how junction IDs, strand IDs, and strand endpoint
crosslinks are identified from that system's files.

## Files

- `junction_strand_network.py`: builds the NetworkX graph, writes a CSV edge
  table, and reports connectivity/functionality checks.
- `plot_junction_strand_network.py`: plots a CSV/graph with junction IDs on
  nodes and strand IDs on edges.
- `src/system_constants.py`: PVA-GLU constants used by the current parser.
- `requirements.txt`: minimal Python dependencies.

## Install

Create or activate a Python environment with:

```bash
pip install -r requirements.txt
```

The required packages are:

- `networkx`
- `numpy`
- `matplotlib`

## Build A CSV

Example for the current PVA-GLU hydrogel parser:

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

The output CSV is an edge table. Each row is one strand connecting two junction
nodes.

Current PVA-GLU CSV columns:

- `edge_id`: row-level strand edge ID
- `source`: source junction ID
- `target`: target junction ID
- `pva_chain_id`: strand/PVA chain ID
- `pva_atom_u`: first PVA endpoint atom
- `pva_atom_v`: second PVA endpoint atom
- `glu_atom_u`: GLU atom bonded to the first PVA endpoint
- `glu_atom_v`: GLU atom bonded to the second PVA endpoint
- `is_self_loop`: `1` if the strand connects back to the same junction

For a fully general version, the chemistry-specific column names can be renamed
to generic labels such as `strand_id`, `strand_endpoint_u`,
`strand_endpoint_v`, `junction_atom_u`, and `junction_atom_v`.

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

The plot uses an abstract graph layout. Node labels are junction IDs and edge
labels are strand IDs, so the figure can be cross-checked directly against the
CSV.

## What To Check

For an ideal loop-free `A2 + B4` network:

- every `B4` junction has degree 4
- every `A2` strand connects two distinct junctions
- there are no self-loop strands
- there are no dangling strands
- there are no unexpected disconnected components
- parallel strand edges are visible in the CSV and graph summary

For the PVA-GLU hydrogel example, the corresponding statement is:

- every GLU junction should have four PVA branches
- every PVA strand should connect two distinct GLU junctions
- the network should be connected
- self-loop and unexpected parallel strand defects should be absent unless they
  are intentionally allowed
