import csv
import os
import gridify_python as gridify
from codetiming import Timer
from concurrent.futures import ProcessPoolExecutor

opt = "FAIL"
opt = 0

dude_dir = "./dude"
out_dir = "run3_mt"

opt = "run3"
save_grid_pdbs = False

os.makedirs(out_dir, exist_ok=True)

def write_pdb(filename, grid):
	idx = 0

	with open(filename, 'w') as file:
		for pt in grid:
			x = pt.x()
			y = pt.y()
			z = pt.z()
			print("HETATM{:>5}  C   PTH {:>5}{:>12.3f}{:>8.3f}{:>8.3f}  0.00  0.00".format(idx, idx, x, y, z), file=file)
			idx += 1

def get_all_combos():
	spacing = 1
	dense_packing = False
	rm_lc_tangent_weight = 1
	rm_lc_proximity_weight = 0.5
	pca_align = False

	if opt == "run1":
		lp_distance_opts = [3, 3.5, 4, 4.5, 5, 5.5]
		lp_ignore_radii_opts = [True]
		scale_radius_opts = [1, 1.5]
		rm_atom_overlaps_opts = [True]
		largest_cluster_only_opts = [False, True]
		rm_lc_cutoff_opts = [0, 1, 2, 3, 4, 5, 6, 7]

	elif opt == "run2":
		lp_distance_opts = [0, 0.5 ,1 ,1.5 ,2 ,2.5 ,3]
		lp_ignore_radii_opts = [False]
		scale_radius_opts = [1, 1.5]
		rm_atom_overlaps_opts = [True]
		largest_cluster_only_opts = [False, True]
		rm_lc_cutoff_opts = [0, 1, 2, 3, 4, 5, 6, 7]

	elif opt == "run3":
		spacing = 0.5
		lp_distance_opts = [0, 0.5 ,1 ,1.5 ,2 ,2.5 ,3]
		lp_ignore_radii_opts = [False]
		scale_radius_opts = [1, 1.5]
		rm_atom_overlaps_opts = [True]
		largest_cluster_only_opts = [False, True]
		rm_lc_cutoff_opts = [0, 1, 2, 3, 4, 5, 6, 7]

	elif opt == 0:
		lp_distance_opts = [2.5, 3, 3.5, 4]
		lp_ignore_radii_opts = [True]
		scale_radius_opts = [1, 1.5]
		rm_atom_overlaps_opts = [True]
		largest_cluster_only_opts = [False, True]
		rm_lc_cutoff_opts = [0, 2, 4, 6, 7, 7.5]	
	else:
		lp_distance_opts = [1, 2, 3]
		lp_ignore_radii_opts = [True, False]
		scale_radius_opts = [1, 1.5]
		rm_atom_overlaps_opts = [True, False]
		largest_cluster_only_opts = [True, False]
		rm_lc_cutoff_opts = [0, 1, 3, 6]


	for lp_distance in lp_distance_opts:
		for lp_ignore_radii in lp_ignore_radii_opts:
			for scale_radius in scale_radius_opts:
				for rm_atom_overlaps in rm_atom_overlaps_opts:
					for largest_cluster_only in largest_cluster_only_opts:
						for rm_lc_cutoff in rm_lc_cutoff_opts:
							yield {
								"lp_distance": lp_distance,
								"lp_ignore_radii": lp_ignore_radii,
								"scale_radius": scale_radius,
								"pca_align": pca_align,
								"spacing": spacing,
								"rm_atom_overlaps": rm_atom_overlaps,
								"largest_cluster_only": largest_cluster_only,
								"dense_packing": dense_packing,
								"rm_lc_cutoff": rm_lc_cutoff,
								"point_radius": spacing / 2,
								"rm_lc_tangent_weight": rm_lc_tangent_weight,
								"rm_lc_proximity_weight": rm_lc_proximity_weight,
							}

def get_stats(protein_chain, ligand_atoms, lp_distance, lp_ignore_radii, rm_atom_overlaps, largest_cluster_only, dense_packing, rm_lc_cutoff, scale_radius, spacing, point_radius, rm_lc_tangent_weight, rm_lc_proximity_weight, pca_align):
	residues = gridify.get_protein_residues_near_ligand(protein_chain, ligand_atoms, lp_distance, lp_ignore_radii)
	site = gridify.get_binding_site(protein_chain, residues)
	grid = gridify.gen_site_grid(protein_chain, site, rm_atom_overlaps, largest_cluster_only, dense_packing, rm_lc_cutoff, scale_radius, spacing, point_radius, rm_lc_tangent_weight, rm_lc_proximity_weight)
	if pca_align:
		grid = gridify.pca_aligned_points(grid)

	ligand = gridify.gen_ligand_geometry(ligand_atoms, 1, pca_align)
	return gridify.calc_grid_ligand_stats(grid, ligand, point_radius, dense_packing), grid

def parse_pdb_and_get_stats(prot_name, ligand_chain, ligand_resid):
	pdb_file = os.path.join(dude_dir, f"{prot_name}.pdb")
	out_file = os.path.join(out_dir, f"{prot_name}.csv")
	print(f"processing protein {pdb_file} with ligand {ligand_resid}-{ligand_chain}")
	with open(out_file, 'w', newline='') as out:
		writer = csv.writer(out)
		writer.writerow(["config idx", "site volume", "ligand volume", "intersection volume", "union volume"])

		frame = gridify.parse_pdb(pdb_file)[0]
		protein = [atom for atom in frame.atoms if atom.kind == "ATOM"]
		protein_chain = [atom for atom in protein if atom.chain == ligand_chain]
		ligand_atoms = [atom for atom in frame.atoms if atom.kind == "HETATM" and atom.chain == ligand_chain and atom.residue == ligand_resid]
		if len(ligand_atoms) == 0:
			return
			#print(f" -- NUM PROTEIN CHAIN ATOMS: {len(protein_chain)}")
			#print(f" -- NUM LIGAND ATOMS: {len(ligand_atoms)}")
		for idx, args in enumerate(get_all_combos()):
			stats, grid = get_stats(protein_chain, ligand_atoms, **args)
			if save_grid_pdbs:
				base_path = os.path.join(out_dir, "grid_pdbs", prot_name)
				os.makedirs(base_path, exist_ok=True)
				write_pdb(os.path.join(base_path, f"option_{idx}.pdb"), grid)
			writer.writerow([idx, stats.site_volume, stats.ligand_volume, stats.intersection_volume, stats.union_volume])

def wrap(row):
        if len(row) == 0:
                return
        prot_name, chain, resid = row
        parse_pdb_and_get_stats(prot_name, chain, resid)

def main():
	i = input("Program init... Click [ENTER] to continue")
	with open(os.path.join(out_dir, "options.txt"), "w") as f:
		for idx, cmds in enumerate(get_all_combos()):
			print(f"{idx}: ", cmds, file=f)

	with open("dude/dude.csv") as csvfile:
		reader = csv.reader(csvfile)
		ProcessPoolExecutor().map(wrap, reader)

if __name__ == '__main__':
    main()