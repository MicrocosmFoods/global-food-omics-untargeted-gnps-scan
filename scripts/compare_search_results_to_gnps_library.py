#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This script takes a list of metabolites of interest (MOI) and...

# Finds if they are in the GNPS library
# Calculates their monoisotopic mass
# Calculates the monoisotopic mass of common adducts
# Pulls GNPS library spectra of the MOI (from downloaded, full GNPS library)
# Searches local mzXML files for all MS2 spectra with precursors within ppm tolerance of MOI +/- adducts
# Calculates similarity scores for the MS2 spectra against the GNPS library spectra
# Outputs a dataframe with the results

import ast
import pandas as pd
from pathlib import Path
from matchms import Spectrum
from matchms.similarity import ModifiedCosine
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
import multiprocessing as mp
from tqdm import tqdm
from matchms.importing import load_from_mgf
from matchms.filtering import remove_peaks_outside_top_k
import pandas as pd
import ast
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import pickle


def load_gnps_library(file_mgf):
    pickle_path = Path(file_mgf).resolve().with_suffix(".pkl")
    if pickle_path.exists():
        with open(pickle_path, "rb") as f:
            return pickle.load(f)
    else:
        print(f"Loading GNPS library from {file_mgf}; fairly slow step...")
        mgf = {i.metadata["spectrum_id"]: i for i in load_from_mgf(file_mgf)}
        with open(pickle_path, "wb") as f:
            pickle.dump(mgf, f)


def add_matchms_spectrum_col(row, intensity="i_tic_norm"):
    """
    Create a matchms Spectrum object from a row of the dataframe."""
    match intensity:
        case "i_tic_norm":
            intensities = row["i_tic_norm"]
        case "i_norm":
            intensities = row["i_norm"]
        case _:
            raise ValueError("Invalid intensity type. Choose 'i_tic_norm' or 'i_norm'.")
    return Spectrum(
        mz=row["mz"],
        intensities=intensities,
        metadata={
            "filename": row["filename"],
            "scan": row["scan"],
            "name": f"{row['filename']}_{row['scan']}",
            "species": row["species"],
            "precursor_mz": row["precmz"],
            "id": f"{row['filename']}_{row['scan']}",
        },
    )


def merge_search_results_and_gnps_clusterinfo(
    molecule_of_interest,
    search_results_path,
    metabolites_of_interest_gnps_library_lookup_path,
    gnps2_results_dir,
):
    """
    Reads the data for a specific molecule of interest from the custom search results directory and
    merges it with cluster information from the GNPS2 results directory.
    Args:
        molecule_of_interest (str): The name of the molecule of interest.
        search_results_path (str): Path to the custom search results directory.
        metabolites_of_interest_gnps_library_lookup_path (str): Path to the metabolites of interest TSV file.
        gnps2_results_dir (str): Path to the GNPS2 results directory.
    """
    df = pd.read_csv(Path(metabolites_of_interest_gnps_library_lookup_path), sep="\t")
    pubchem_cid = df[df["metabolite_of_interest"] == molecule_of_interest][
        "pubchem_cid"
    ].values[0]
    df = Path(
        search_results_path,
        f"{pubchem_cid}_{molecule_of_interest}_mslevel-2.feather",
    ).resolve()
    df = pd.read_feather(df)
    df["filename_standardized"] = (
        df["filename"].str.split("/").str[-1].str.removesuffix(".mzXML")
    )
    clusterinfo_df = pd.read_csv(
        Path(
            gnps2_results_dir,
            "nf_output/clustering/clusterinfo.tsv",
        ),
        sep="\t",
    )
    clusterinfo_df["filename_standardized"] = (
        clusterinfo_df["#Filename"].str.split("/").str[-1].str.removesuffix(".mzML")
    )
    clusterinfo_df = clusterinfo_df[["filename_standardized", "#Scan", "#ClusterIdx"]]
    # make scan an int
    clusterinfo_df["#Scan"] = clusterinfo_df["#Scan"].astype(int)
    # change column names
    clusterinfo_df.columns = ["filename_standardized", "scan", "gnps_cluster"]
    df = pd.merge(df, clusterinfo_df, on=["filename_standardized", "scan"], how="left")
    df["sample_spectrum"] = df.apply(
        lambda x: add_matchms_spectrum_col(x, intensity="i_tic_norm"), axis=1
    )
    df.drop(columns=["i_norm", "i_tic_norm", "mz"], inplace=True)
    return df


def compare_to_gnps_lib_spectrum(
    df,
    ccmslib_spectrum,
    ccmslib_metadata,
    algorithm="CosineGreedy",
    score_threshold=0.5,
    k_top_peaks=20,
    mz_window=5,
):
    """
    Note: k_top_peaks "Removes all peaks which are not within mz_window of at least one of the k highest intensity peaks of the spectrum.
    see: https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html#matchms.filtering.remove_peaks_outside_top_k
    """
    match algorithm:
        case "CosineGreedy":
            similarity_measure = CosineGreedy()
        case "ModifiedCosine":
            similarity_measure = ModifiedCosine()
        case _:
            raise ValueError(
                "Invalid algorithm. Choose 'CosineGreedy' or 'ModifiedCosine'."
            )
    new_df = df.copy()
    filtered_sample_spectra = (
        new_df["sample_spectrum"]
        .apply(
            lambda x: remove_peaks_outside_top_k(
                spectrum_in=x, k=k_top_peaks, mz_window=mz_window
            )
        )
        .to_list()
    )
    filtered_ccmslib_spectrum = remove_peaks_outside_top_k(
        spectrum_in=ccmslib_spectrum, k=k_top_peaks, mz_window=mz_window
    )
    scores = calculate_scores(
        [filtered_ccmslib_spectrum], filtered_sample_spectra, similarity_measure
    )
    array_data = scores.to_array().flatten()
    df_scores = pd.DataFrame(array_data)
    df_with_scores = new_df.merge(df_scores, left_index=True, right_index=True)
    new_df = df_with_scores[
        (df_with_scores[f"{algorithm}_score"] > score_threshold)
        & (df_with_scores["gnps_cluster"].notna())
    ]
    new_df = new_df.sort_values(
        [f"{algorithm}_score", f"{algorithm}_matches"], ascending=False
    )
    new_df["ccmslib_id"] = ccmslib_spectrum.metadata["spectrum_id"]
    new_df["ccmslib_spectrum"] = ccmslib_spectrum
    # ccmslib_metadata = {'spectrum_id': 'CCMSLIB00005463716', 'source_file': 'f.ericf74/reference_spectra/reference_spectra/49ada27d-7c4f-4018-bdf5-08b17b587d78.mgf;', 'task': '2da3086a339648e3b964d4499ab86770', 'scan': '1', 'ms_level': '2', 'library_membership': 'GNPS-LIBRARY', 'spectrum_status': '1', 'peaks_json': 'null', 'splash': 'null-null-null-null', 'submit_user': 'ericf74', 'Compound_Name': 'Tryptophane-emf', 'Ion_Source': 'LC-ESI', 'Compound_Source': 'Crude', 'Instrument': 'Orbitrap', 'PI': 'emf', 'Data_Collector': 'N/A', 'Adduct': 'M+H', 'Scan': '-1', 'Precursor_MZ': '205.097', 'ExactMass': '204.09', 'Charge': '1', 'CAS_Number': 'N/A', 'Pubmed_ID': ' ', 'Smiles': 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N', 'INCHI': 'InChI=1S/C11H12N2O2/c12-9(11(14)15)5-7-6-13-10-4-2-1-3-8(7)10/h1-4,6,9,13H,5,12H2,(H,14,15)/t9-/m0/s1', 'INCHI_AUX': 'QIVBCDIJIAJPQS-VIFPVBQESA-N', 'Library_Class': '3', 'SpectrumID': 'CCMSLIB00005463716', 'Ion_Mode': 'Positive', 'create_time': '2021-03-18 20:16:09.0', 'task_id': '90fc233158e44408bca351458ccfebb6', 'user_id': 'null'}
    # all at once add meta columns for Library_Class, Compound_Source, Instrument, spectrum_status, library_membership,
    new_df["library_class"] = ccmslib_metadata["Library_Class"]
    new_df["compound_source"] = ccmslib_metadata["Compound_Source"]
    new_df["instrument"] = ccmslib_metadata["Instrument"]
    new_df["spectrum_status"] = ccmslib_metadata["spectrum_status"]
    new_df["library_membership"] = ccmslib_metadata["library_membership"]
    return new_df


def process_row(
    row,
    search_results_path,
    moi_lookup_path,
    gnps2_results_dir,
    output_path,
    algorithm,
    score_threshold,
    k_top_peaks,
    mz_window,
):
    moi_results = []
    # Generate search results
    df = merge_search_results_and_gnps_clusterinfo(
        molecule_of_interest=row["metabolite_of_interest"],
        search_results_path=search_results_path,
        metabolites_of_interest_gnps_library_lookup_path=moi_lookup_path,
        gnps2_results_dir=gnps2_results_dir,
    )
    for spectrum_id, data in row["submatch_positive"].items():
        try:
            df2 = compare_to_gnps_lib_spectrum(
                df=df,
                ccmslib_spectrum=data["spectrum"],
                ccmslib_metadata=data["metadata"],
                algorithm=algorithm,
                score_threshold=score_threshold,
                k_top_peaks=k_top_peaks,
                mz_window=mz_window,
            )
            moi_results.append(df2)
        except Exception as e:
            print(f"Error processing spectrum {spectrum_id}: {e}")
            continue
    dfend = pd.concat(moi_results) if moi_results else None
    # Pickle write dfend
    with open(output_path / f"{row['metabolite_of_interest']}.pkl", "wb") as f:
        pickle.dump(dfend, f)


def run_in_parallel(
    gnps_mgf,
    metabolites_of_interest_gnps_library_lookup_path,
    search_results_path,
    moi_lookup_path,
    gnps2_results_dir,
    output_path,
    algorithm,
    score_threshold,
    k_top_peaks,
    mz_window,
):
    # Read data
    ccmslib = load_gnps_library(gnps_mgf)
    moi_df = pd.read_csv(metabolites_of_interest_gnps_library_lookup_path, sep="\t")
    moi_df = moi_df[["metabolite_of_interest", "submatch_positive"]]
    moi_df["submatch_positive"] = moi_df["submatch_positive"].apply(ast.literal_eval)
    # Convert `submatch_positive` column to a dictionary
    moi_df["submatch_positive"] = moi_df["submatch_positive"].apply(
        lambda x: {
            i["spectrum_id"]: {"spectrum": ccmslib[i["spectrum_id"]], "metadata": i}
            for i in x
            if i["spectrum_id"] in ccmslib
        }
    )
    del ccmslib
    n_workers = min(cpu_count(), len(moi_df))  # Number of workers
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                process_row,
                row,
                search_results_path,
                moi_lookup_path,
                gnps2_results_dir,
                output_path,
                algorithm,
                score_threshold,
                k_top_peaks,
                mz_window,
            ): idx
            for idx, row in moi_df.iterrows()
        }
        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Processing rows"
        ):
            try:
                future.result()
            except Exception as e:
                print(f"Error processing row {futures[future]}: {e}")
