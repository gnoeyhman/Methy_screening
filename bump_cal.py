import logging
from Bio import PDB
from pymol import cmd
import os
import glob

# Set up logging
logging.basicConfig(filename='execution.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def renum_h(ref_molecule):
    cmd.load(indir + '/' + ref_molecule, "ref_mol")
    cmd.remove("hydrogens")
    cmd.h_add("all")
    cmd.save(indir + '/' + f"{ref_molecule[:-4]}_h.pdb", "all")
    cmd.delete("all")
    cmd.reinitialize

def convert_hydrogens_to_carbons(ref_molecule, indir):
    ref_mol_renum = ref_molecule[:-4] + '_h.pdb'
    parser = PDB.PDBParser()
    structure = parser.get_structure('structure', indir + '/' + ref_mol_renum)

    for model in structure:
        for chain in model:
            for residue in chain:
                hydrogen_atoms = [atom for atom in residue if atom.element == 'H']
                for i, hydrogen_atom in enumerate(hydrogen_atoms, start=1):
                    new_atom = hydrogen_atom.copy()
                    new_atom.element = 'C'

                    # Detach the original hydrogen atom
                    residue.detach_child(hydrogen_atom.id)

                    # Add the new carbon atom
                    residue.add(new_atom)

                    # Write the modified structure to a new PDB file
                    output_pdb_file = os.path.join(indir, f"converted_{i}.pdb")
                    io = PDB.PDBIO()
                    io.set_structure(structure)
                    io.save(output_pdb_file)

                    # Revert the change to the residue for the next conversion
                    residue.detach_child(new_atom.id)
                    residue.add(hydrogen_atom)

'''
http://pymolwiki.org/index.php/show_bumps

(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

def show_bumps(selection='(all)', name='bump_check', quiet=1):
    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
    for state in range(1, 1 + cmd.count_states('%' + name)):
        cmd.sculpt_activate(name, state)
        strain = cmd.sculpt_iterate(name, state, cycles=0)
        if not int(quiet):
            print('VDW Strain in state %d: %f' % (state, strain))
    cmd.show_as('cgo', name)

def convert_0():
    cmd.load(indir + '/' + apo_str, "apo")
    cmd.load(indir + '/' + ref_molecule, "ref_mol")

    selection = "apo or ref_mol"
    name = "test_0"
    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)
    cmd.sculpt_activate(name, 1)
    strain = cmd.sculpt_iterate(name, 1, cycles=0)
    if not int(1):
        print('VDW Strain in state %d: %f' % (1, strain))
    cmd.show_as('cgo', name)
    print(f"vdw strain in reference: {strain}")
    cmd.delete("all")
    cmd.reinitialize

def count_files_with_prefix(indir):
    # Initialize count to 0
    count = 0
    # Iterate over files in the specified directory
    for filename in os.listdir(indir):
        # Check if the file starts with the specified prefix
        if filename.startswith("converted_"):
            # Increment count if the condition is met
            count += 1
            print(count)
    return count
    print(count)

def calculate_bump_net(indir, outdir):
    count = 0
    # Iterate over files in the specified directory
    for filename in os.listdir(indir):
        # Check if the file starts with the specified prefix
        if filename.startswith("converted_"):
            # Increment count if the condition is met
            count += 1

    ref_mol_renum = ref_molecule[:-4] + '_h.pdb'

    for i in range(1, count + 1):
        cmd.load(indir + '/' + apo_str, "apo")
        cmd.load(indir + '/' + ref_mol_renum, "ref_mol")
        cmd.load(indir + '/' + f"converted_{i}.pdb", f"test_{i}")
        cmd.create(f"mol_vis_{i}", f"test_{i}")
        cmd.create("apo_vis", "apo")
        cmd.remove(f"hydrogens and (mol_vis_{i} or apo_vis)")

        select_ref = f"apo or ref_mol"
        selection = f"apo or test_{i}"
        select_vis = f"apo_vis or mol_vis_{i}"
        name_ref = "test_0"
        name = f"bump_{i}"
        name_vis = f"vis_{i}"
        cmd.delete(name_ref)
        cmd.delete(name)
        cmd.delete(name_vis)
        cmd.create(name_ref, select_ref, zoom=0)
        cmd.create(name, selection, zoom=0)
        cmd.create(name_vis, select_vis, zoom=0)

        cmd.set('sculpt_vdw_vis_mode', 1, name_ref)
        cmd.set('sculpt_vdw_vis_mode', 1, name)
        cmd.set('sculpt_vdw_vis_mode', 1, name_vis)
        cmd.set('sculpt_field_mask', 0x020)
        with open(outdir + '/' + log_file, 'a') as f:
            cmd.sculpt_activate(name_ref, 1)
            cmd.sculpt_activate(name, 1)
            cmd.sculpt_activate(name_vis, 1)
            strain_ref = cmd.sculpt_iterate(name_ref, 1, cycles=0)
            strain = cmd.sculpt_iterate(name, 1, cycles=0)
            strain_vis = cmd.sculpt_iterate(name_vis, 1, cycles=0)
            cmd.show_as('cgo', name_ref)
            cmd.show_as('cgo', name)
            cmd.show_as('cgo', name_vis)
            #print(f"vdw strain in test_{i}: {strain}")
            print(f"net vdw strain in test_{i}: {strain - strain_ref:.3f}")
            f.write(f"net vdw strain in test_{i}: {strain - strain_ref:.3f}\n")
            cmd.save(outdir + '/' + f"bump_{i}.pse")
            cmd.delete("all")
            cmd.reinitialize

def sort_file(log_file, cut_off):
    with open(outdir + '/' + log_file, 'r') as f:
        lines = f.readlines()

    with open(outdir + '/' + log_file, 'a') as f:
        f.write("\n---Possible candidate---\n")
        for line in lines:
            parts = line.split(':')
            if len(parts) == 2:
                bump_value = float(parts[1].strip())
                if bump_value < float(vdw_cut_off):
                    sorts = line.split(' ')
                    sort_name = sorts[4].strip()
                    hydrogen_sort = sort_name.split('_')
                    hydrogen_number = hydrogen_sort[-1]
                    f.write('H' + hydrogen_number[:-1] + '/')
                else:
                    continue

def visulaize(indir, outdir, apo_str, ref_molecule, select_cut_off):
    ref_mol_renum = ref_molecule[:-4] + '_h.pdb'

    cmd.load(indir + '/' + apo_str, "apo")
    cmd.load(indir + '/' + ref_mol_renum, "ref_mol")
    with open(outdir + '/' + log_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if '/' in line:
                h_sort = line.split('/')[:-1]
                print(h_sort)
                for i in h_sort:
                    cmd.select(i, f"apo and (br. all within {select_cut_off} of ref_mol and name {i})")
                    cmd.show("sticks", i)
                    cmd.color("red", i, quiet=1)
        cmd.center("ref_mol")
        cmd.save(outdir + '/' + "visual.pse")

def main():
    os.makedirs(outdir, exist_ok=True)
    logging.info("Execution started.")

    # convert reference as convert_0
    renum_h(ref_molecule)
    convert_hydrogens_to_carbons(ref_molecule, indir)

    # convert reference as convert_0
    convert_0()

    # calculate bump
    calculate_bump_net(indir, outdir)
    sort_file(log_file, vdw_cut_off)
    visulaize(indir, outdir, apo_str, ref_molecule, select_cut_off)

    logging.info("Execution completed.")

if __name__ == "__main__":
    indir = "/Users/gnoeyhman/computation/DeGrado_Lab/carbon_screening/03_04_2024_final/input_2"
    outdir = "/Users/gnoeyhman/computation/DeGrado_Lab/carbon_screening/03_04_2024_final/output"

    log_file = "log_bump.txt"
    vdw_cut_off = 0.3

    select_cut_off = 3.8

    apo_str = "apo_2.pdb"
    ref_molecule = "bas.pdb"

    main()
