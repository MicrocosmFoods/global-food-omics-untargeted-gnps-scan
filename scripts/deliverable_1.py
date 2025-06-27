# %%
# conda create -n py312 python=3.12
# conda activate py312
# pip install requests pubchempy pandas rdkit ipython ipykernel

import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs, MolFromSmiles
from IPython.display import display, HTML
import base64
from io import BytesIO
import time
import requests
import json
from pyopenms import MSExperiment, MzXMLFile
from pathlib import Path


def get_compound_info(cid, retries=3, delay=0):
    try:
        # Add delay to avoid hitting rate limit
        time.sleep(1)
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


def img_to_base64(mol):
    # Generate image from RDKit molecule
    img = Draw.MolToImage(mol)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return f"data:image/png;base64,{img_str}"


def gendata(gnps_library_json_path, output_path):
    if Path(output_path).exists():
        df = pd.read_pickle(output_path)
        return df
    Path(output_path).parent.mkdir(
        parents=True, exist_ok=True
    )  # Ensure the directory exists
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
    df["getcompoundinfotemp"] = df["pubchem_cid"].apply(
        lambda cid: get_compound_info(cid) if pd.notnull(cid) else None
    )
    # Add columns with chemical structure images, molecular formula, monoisotopic mass, and InChI Key
    df["Chemical Structure Image"] = df["getcompoundinfotemp"].apply(
        lambda info: img_to_base64(Chem.MolFromSmiles(info[0])) if info else None
    )
    df["Molecular Formula"] = df["getcompoundinfotemp"].apply(
        lambda info: info[1] if info else None
    )
    df["Monoisotopic Mass"] = df["getcompoundinfotemp"].apply(
        lambda info: round(float(info[2]), 4) if info and info[2] is not None else None
    )
    df["InChI Key"] = df["getcompoundinfotemp"].apply(
        lambda info: info[3] if info else None
    )
    df["rdkit_compound"] = df["getcompoundinfotemp"].apply(
        lambda info: Chem.MolFromSmiles(info[0]) if info else None
    )
    df.drop(columns=["getcompoundinfotemp"], inplace=True)
    with open(gnps_library_json_path, "r") as f:
        gnps_library = json.load(f)
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
            and i["Ion_Mode"] in [" Positive", "Positive", "positive", "Positive-20eV"]
        ]
    )
    df["submatch_negative"] = df["InChI Key"].apply(
        lambda inchi_key: [
            i
            for i in gnps_library
            if i["INCHI_AUX"][:14] == inchi_key[:14]
            and i["Ion_Mode"] in ["Negative", "negative", " Negative"]
        ]
    )
    df = df.sort_values(by="Monoisotopic Mass", ascending=False)
    df.to_pickle(output_path)
    return df


def render_html_table_with_images(df):
    html = """
    <details>
      <summary>Click to show table of Metabolites of Interest</summary>
      <table style='border: 1px solid black; border-collapse: collapse;'>
        <tr>
            <th>Metabolite of Interest</th>
            <th>InChI Key</th>
            <th>PubChem CID</th>
            <th>Formula</th>
            <th>Monoisotopic Mass</th>
            <th>Chemical Structure</th>
            <th>GNPS Library Entries (Positive)</th>
            <th>GNPS Library Entries (Negative)</th>
            <th>All Adducts (Positive)</th>
            <th>All Ion Sources (Positive)</th>
            <th>All Instruments (Positive)</th>
            <th>All Compound Sources (Positive)</th>
        </tr>
    """

    for _, row in df.iterrows():
        html += f"<tr><td>{row['metabolite_of_interest']}</td>"
        html += f"<td>{row['InChI Key']}</td>"
        html += f"<td><a href='https://pubchem.ncbi.nlm.nih.gov/compound/{row['pubchem_cid']}' target='_blank'>{row['pubchem_cid']}</a></td>"
        html += f"<td>{row['Molecular Formula']}</td>"
        html += f"<td>{row['Monoisotopic Mass']}</td>"

        # Chemical structure image
        if row["Chemical Structure Image"]:
            html += (
                f"<td><img src='{row['Chemical Structure Image']}' width='100'></td>"
            )
        else:
            html += "<td>No structure available</td>"

        # GNPS Library Match Counts
        positive_count = len(row.get("submatch_positive", []))
        negative_count = len(row.get("submatch_negative", []))
        html += f"<td>{positive_count}</td>"
        html += f"<td>{negative_count}</td>"

        # Extract unique attributes
        for attribute in ["Adduct", "Ion_Source", "Instrument", "Compound_Source"]:
            unique_values = set(
                match.get(attribute, "") for match in row.get("submatch_positive", [])
            )
            html += f"<td>{', '.join(unique_values)}</td>"

        html += "</tr>"

    html += "</table></details>"
    return html
