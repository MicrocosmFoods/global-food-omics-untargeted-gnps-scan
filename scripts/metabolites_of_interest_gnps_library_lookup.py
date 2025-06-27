import pandas as pd
import pubchempy as pcp
from rdkit import Chem
import time
import requests
import json
from pathlib import Path
import ast


def get_compound_info(cid, retries=3, delay=5):
    try:
        # Add delay to avoid hitting rate limit
        time.sleep(0.2)
        compound = pcp.Compound.from_cid(cid)
        smiles = compound.isomeric_smiles
        formula = compound.molecular_formula
        monoisotopic_mass = compound.monoisotopic_mass
        inchi_key = compound.inchikey
        return smiles, formula, monoisotopic_mass, inchi_key
    except requests.exceptions.RequestException as e:
        if retries > 0:
            print(f"Error fetching CID {cid}. Retrying in {delay} seconds...")
            time.sleep(delay)
            return get_compound_info(cid, retries - 1, delay)
        else:
            print(f"Failed to fetch CID {cid} after multiple attempts: {e}")
            return None, None, None, None


def download_gnps_library(json_path, force=False):
    json_path.parent.mkdir(parents=True, exist_ok=True)
    if not force and json_path.exists():
        print("skipping download of gnpslibrary.json")
        return
    # URL to download from
    url = "https://external.gnps2.org/gnpslibraryjson"
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(json_path, "wb") as f:
            # Iterate over the content in chunks
            for chunk in response.iter_content(chunk_size=8192):
                _ = f.write(chunk)
        print("File downloaded successfully.")
    else:
        print(f"Failed to download the file. Status code: {response.status_code}")


def create_df(json_path, outpath):
    download_gnps_library(json_path)
    metabolites_of_interest_gnps_library_lookup_path = outpath
    if metabolites_of_interest_gnps_library_lookup_path.exists():
        df = pd.read_csv(metabolites_of_interest_gnps_library_lookup_path, sep="\t")
        # ast ensures that the submatch_positive and submatch_negative columns are lists
        df["submatch_positive"] = df["submatch_positive"].apply(ast.literal_eval)
        df["submatch_negative"] = df["submatch_negative"].apply(ast.literal_eval)
        return df
    data = {
        "metabolite_of_interest": [
            "hyodeoxycholic acid",
            "Conjugated linoleic acid",
            "Indole-3-lactic acid",
            "Tryptophan",
            "Indopropionic acid",
            "Indole-3-propionic acid",
            "4-hydroxyphenyllactic acid",
            "Hippuric acid",
            "Serotonin",
            "3-3-PPA",
            "D-phenyllactic acid",
            "Glutamate",
            "alpha-hydroxyisocaproate",
            "2-hydroxy-3-methylvalerate",
            "Succinate",
            "Histamine",
            "3-hydroxybutyric acid",
            "GABA",
            "Lactic acid",
            "TMAO",
            "Propionate",
            "TMA",
            "Acetate",
            "alpha-hydroxyisovalerate",
        ],
        "pubchem_cid": [
            5283820,
            5282796,
            92904,
            6305,
            3744,
            3744,
            9378,
            464,
            5202,
            91,
            444718,
            33032,
            92779,
            10796774,
            160419,
            774,
            441,
            119,
            612,
            1145,
            104745,
            1146,
            175,
            99823,
        ],
    }
    df = pd.DataFrame(data)
    compound_info = df["pubchem_cid"].apply(get_compound_info)
    df[["SMILES", "Molecular Formula", "Monoisotopic Mass", "InChI Key"]] = (
        pd.DataFrame(compound_info.tolist(), index=df.index)
    )
    df.to_csv("metabolites_of_interest_gnps_library_lookup.tsv", sep="\t", index=False)
    with open(json_path, "r") as f:
        gnps_library = json.load(f)
    # The JSON file doesn't have all of the INHCHI keys that the web data does, so where it is missing, we will try to add it via calculating it from the SMILES or InChI
    for i in gnps_library:
        if i["INCHI_AUX"] != "N/A":
            continue
        elif i["Smiles"] != "N/A":
            try:
                mol = Chem.MolFromSmiles(i["Smiles"])
                i["INCHI_AUX"] = Chem.MolToInchiKey(mol)
            except:
                pass
        elif i["INCHI"] != "N/A":
            try:
                mol = Chem.MolFromInchi(i["InChI"])
                i["INCHI_AUX"] = Chem.MolToInchiKey(mol)
            except:
                pass
    df["submatch_positive"] = df["InChI Key"].apply(
        lambda inchi_key: [
            i
            for i in gnps_library
            if i["INCHI_AUX"][:14] == inchi_key[:14]
            and i["Ion_Mode"]
            in ["", "Positive", "positive", "Positive-20eV", " Positive", "N/A"]
        ]
    )
    df["submatch_negative"] = df["InChI Key"].apply(
        lambda inchi_key: [
            i
            for i in gnps_library
            if i["INCHI_AUX"][:14] == inchi_key[:14]
            and i["Ion_Mode"] in ["", "negative", "N/A", "Negative", " Negative"]
        ]
    )
    # ast ensures that the submatch_positive and submatch_negative columns are lists
    df["submatch_positive"] = df["submatch_positive"].apply(str)
    df["submatch_negative"] = df["submatch_negative"].apply(str)
    df.to_csv(outpath, sep="\t", index=False)
