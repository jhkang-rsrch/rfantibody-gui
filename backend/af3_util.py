from Bio import PDB
from pathlib import Path
import json
from Bio.Data import IUPACData

def prepare_af3_json_input(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    input_dir_cifs = list(input_dir.glob("*.cif"))

    for cif_file in input_dir_cifs:
        output_json_file = output_dir / f"{cif_file.stem}.json"
        json_maker(cif_file, output_json_file)


def json_maker(input_cif, output_json):
    input_cif = Path(input_cif)
    output_json = Path(output_json)

    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("prot", input_cif)

    sequences = []
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    resname = residue.get_resname().capitalize()
                    seq += IUPACData.protein_letters_3to1.get(resname, "X")
            sequences.append({
                "protein": {
                    "id": chain.id,
                    "sequence": seq,
                    # "unpairedMsa": None,
                    # "pairedMsa": None,
                    # "templates": [{"mmcifPath": str(input_cif),
                    #                "queryIndices": [i for i in range(len(seq))],
                    #                "templateIndices": [i for i in range(len(seq))]
                    #               }]
                }
            })

    af3_json = {
        "name": input_cif.stem,
        "modelSeeds": [42],
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 4, ## ★ Original was 1 ★
    }

    with open(output_json, "w") as f:
        json.dump(af3_json, f, indent=2)

    print(f"✅ JSON saved to: {output_json}")


def prepare_af3_cif_input(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    parser = PDB.PDBParser()

    input_dir_pdbs = list(input_dir.glob("*.pdb"))
    for idx, pdb_file in enumerate(input_dir_pdbs):
        structure = parser.get_structure(f"X_{idx}", pdb_file)
        io = PDB.MMCIFIO()
        io.set_structure(structure)
        io.save((output_dir / f"input_{idx+1}.cif").as_posix())
