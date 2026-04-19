"""Shared structural constants for the hydrogel system."""

PVA_ATOMS_PER_MONOMER = 10
GLU_ATOMS_PER_MOLECULE = 31
# Chosen to keep total atom count < 80,000 at n=25 with n_pva = 2 * n_glu.
# Total atoms = n_glu * (2 * PVA_ATOMS_PER_MONOMER * 25 + GLU_ATOMS_PER_MOLECULE)
# => 150 * 531 = 79,650 atoms.
DEFAULT_N_GLU = 150
