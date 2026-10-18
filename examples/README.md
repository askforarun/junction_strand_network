# Production topology examples

Each directory contains one representative fully crosslinked production
network exported from the corrected fully hydrolyzed PVA--GLU workflow.
These CSVs and reports are self-contained for graph analysis and plotting;
the large atomic LAMMPS coordinate files are not bundled here.

| Directory | Strand length | Junctions | Strands | Source configuration ID |
|---|---:|---:|---:|---|
| `production_n17` | 17 | 1350 | 2700 | `1180c166f2e323e3e0efb9efd5931fbe` |
| `production_n21` | 21 | 1050 | 2100 | `90a60f1bf1d090ebe85d2fc4030c00c7` |
| `production_n25` | 25 | 900 | 1800 | `4e4d5e4ec34deb6ef6af84bbbca1f5e4` |

Each export used the source `data.crosslinked_updated` with its corresponding
`PVA.txt` and `GLU.txt`, and reactive atom types 9 and 10. Configuration IDs
document provenance; they are not runtime dependencies or paths.

Files in each directory:

- `junction_edges.csv`: one row per strand, including both endpoint bond pairs.
- `PVA.txt`: two reactive atom IDs per strand.
- `GLU.txt`: four reactive atom IDs per junction.
- `junction_strand_network_report.txt`: connectivity and graph diagnostics.

From the repository root:

```bash
python plot_junction_strand_network.py \
  --csv examples/production_n25/junction_edges.csv \
  --out production_n25.png --layout spring
```

The root-level atomic example is a separate, smaller demonstration with
150 junctions and 300 strands. Do not mix its atomic data with the production
index files. A graph diagram shows connectivity, not a periodic molecular
snapshot, residual strain, or geometric spatial percolation.
