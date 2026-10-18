"""Synthetic full-conversion A2+B4 networks and standalone CLI checks."""
import contextlib
import io
import itertools
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from junction_strand_network import (
    build_junction_strand_network, build_pva_edges_and_chain_map,
    summarize_junction_strand_network,
)
from export_and_plot_junction_strand_network import main
from plot_junction_strand_network import read_junction_strand_network_data


class NetworkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.data = self.root / "data.crosslinked"
        self.pva, self.glu = self.root / "PVA.txt", self.root / "GLU.txt"
        self.pva.write_text("\n".join(map(str, range(1, 21))))
        self.glu.write_text("\n".join(map(str, range(21, 41))))
        sites = {node: iter(range(21+4*node, 25+4*node)) for node in range(5)}
        self.bonds = []
        for chain, (u, v) in enumerate(itertools.combinations(range(5), 2)):
            self.bonds.extend([(2*chain+1, next(sites[u])), (2*chain+2, next(sites[v]))])
        self.write_data(self.bonds)

    def write_data(self, bonds):
        atoms = "\n".join(f"{i} 1 {9 if i<=20 else 10} 0 0 0 0" for i in range(1,41))
        lines = "\n".join(f"{i} 18 {u} {v}" for i,(u,v) in enumerate(bonds,1))
        self.data.write_text(f"Synthetic network\n\n40 atoms\n{len(bonds)} bonds\n\nAtoms # full\n\n{atoms}\n\nBonds\n\n{lines}\n")

    def build(self):
        return build_junction_strand_network(str(self.data),str(self.pva),str(self.glu),9,10)

    def test_complete_network_and_cli(self):
        summary = summarize_junction_strand_network(self.build())
        self.assertEqual(summary["n_junction_nodes"],5)
        self.assertEqual(summary["n_strand_edges"],10)
        self.assertTrue(summary["has_4_graph_connectivity"])
        output = self.root / "output"
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(["--data",str(self.data),"--pva",str(self.pva),
                                   "--glu",str(self.glu),"--output-dir",str(output),"--skip-plot"]),0)
        self.assertEqual(read_junction_strand_network_data(str(output/"junction_edges.csv")).number_of_edges(),10)
        self.assertIn("Connected: YES",(output/"junction_strand_network_report.txt").read_text())

    def test_missing_whole_strand_rejected(self):
        self.write_data(self.bonds[2:])
        with self.assertRaisesRegex(ValueError,"no formed crosslinks"):
            self.build()

    def test_reused_junction_site_rejected(self):
        bonds = list(self.bonds)
        bonds[2] = (bonds[2][0],bonds[0][1])
        self.write_data(bonds)
        with self.assertRaisesRegex(ValueError,"GLU reactive sites"):
            self.build()

    def test_duplicate_and_fractional_ids_rejected(self):
        for content in ("1\n1\n", "1\n2.5\n"):
            self.pva.write_text(content)
            with self.assertRaises(ValueError):
                build_pva_edges_and_chain_map(str(self.pva))

    def test_bundled_corrected_chemistry_example(self):
        graph = build_junction_strand_network(str(ROOT/"data.crosslinked"),
                    str(ROOT/"PVA.txt"),str(ROOT/"GLU.txt"),9,10)
        self.assertEqual(graph.number_of_nodes(),150)
        self.assertEqual(graph.number_of_edges(),300)

    def test_production_csvs_and_reactive_ids(self):
        for n, nodes, edges in ((17,1350,2700),(21,1050,2100),(25,900,1800)):
            folder = ROOT / "examples" / f"production_n{n}"
            with self.subTest(n=n):
                graph = read_junction_strand_network_data(str(folder/"junction_edges.csv"))
                self.assertEqual(graph.number_of_nodes(),nodes)
                self.assertEqual(graph.number_of_edges(),edges)
                self.assertTrue(all(degree==4 for _,degree in graph.degree()))
                endpoints, sites = [], []
                for _,_,data in graph.edges(data=True):
                    endpoints.extend([data["pva_atom_u"],data["pva_atom_v"]])
                    sites.extend([data["glu_atom_u"],data["glu_atom_v"]])
                self.assertEqual(sorted(endpoints),sorted(map(int,(folder/"PVA.txt").read_text().split())))
                self.assertEqual(sorted(sites),sorted(map(int,(folder/"GLU.txt").read_text().split())))


if __name__ == "__main__":
    unittest.main()
