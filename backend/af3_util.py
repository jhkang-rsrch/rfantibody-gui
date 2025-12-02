from Bio import PDB
from pathlib import Path
import json
from Bio.Data import IUPACData


STAGE_SUFFIX = {
    "rfdiffusion": "_rfd",      # RFdiffusion 출력
    "mpnn": "_mpnn",            # ProteinMPNN 출력  
    "cif": "_cif",              # CIF 변환 출력
    "af3input": "_af3in",       # AF3 JSON 입력
}


def prepare_af3_json_input(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    input_dir_cifs = sorted(input_dir.glob("*.cif"))  # 정렬하여 일관된 순서 보장

    for cif_file in input_dir_cifs:
        # 단계 접미사 추가: ab_des_0_mpnn_cif.cif -> ab_des_0_mpnn_cif_af3in.json
        output_json_file = output_dir / f"{cif_file.stem}{STAGE_SUFFIX['af3input']}.json"
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
        "version": 3, ## ★ Original was 1 ★
    }

    with open(output_json, "w") as f:
        json.dump(af3_json, f, indent=2)

    print(f"✅ JSON saved to: {output_json}")


def prepare_af3_cif_input(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    parser = PDB.PDBParser(QUIET=True)

    input_dir_pdbs = sorted(input_dir.glob("*.pdb"))  # 정렬하여 일관된 순서 보장
    for idx, pdb_file in enumerate(input_dir_pdbs):
        structure = parser.get_structure(pdb_file.stem, pdb_file)
        io = PDB.MMCIFIO()
        io.set_structure(structure)
        # 단계 접미사 추가: ab_des_0_mpnn.pdb -> ab_des_0_mpnn_cif.cif
        output_cif = output_dir / f"{pdb_file.stem}{STAGE_SUFFIX['cif']}.cif"
        io.save(output_cif.as_posix())
        print(f"✅ CIF saved: {pdb_file.name} -> {output_cif.name}")
