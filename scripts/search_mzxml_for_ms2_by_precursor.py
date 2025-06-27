import pandas as pd
from pathlib import Path
from functools import partial
from tqdm import tqdm
import multiprocessing
from pyteomics import mass
from massql import msql_fileloading
import argparse

parser = argparse.ArgumentParser(description="Process mzXML files and metabolites data")
parser.add_argument(
    "--basedir", type=str, required=True, help="Base directory for input files"
)
parser.add_argument(
    "--subdir", type=str, required=True, help="Subdirectory under base directory"
)
parser.add_argument(
    "--ppm-tolerance", type=float, default=3.0, help="PPM tolerance for mass matching"
)
parser.add_argument(
    "--cpus",
    type=int,
    default=20,
    help="Number of worker processes for parallel processing",
)
parser.add_argument("--level", type=int, default=2, help="Level of MS data to search")
parser.add_argument(
    "--output-dir",
    type=str,
    default="massql_search/search_results",
    help="Output directory for search results",
)
parser.add_argument(
    "--metabolites-tsv", type=str, required=True, help="Path to metabolites data"
)
parser.add_argument(
    "--adducts-tsv", type=str, required=True, help="Path to mass adducts data"
)


def calc_expected_mass(molecular_mass, species_mass, species_mult):
    return (molecular_mass * species_mult) + species_mass


def calc_ppm(mass, ppm_tolerance):
    return mass * ppm_tolerance / 1e6


def process_file(
    file, metabolites, ion_mass_map, ppm_tolerance, basedir, level, output_dir
):
    ms1_df, ms2_df = msql_fileloading.load_data(str(file))
    relative_filename = file.relative_to(basedir)
    for _, row in metabolites.iterrows():
        for species, (species_mass, species_mult) in ion_mass_map.items():
            expected_mass = calc_expected_mass(
                row["monoisotopic_mass"], species_mass, species_mult
            )
            lower = expected_mass - calc_ppm(expected_mass, ppm_tolerance)
            upper = expected_mass + calc_ppm(expected_mass, ppm_tolerance)
            if level == 1:
                results_df = ms1_df[(ms1_df["mz"] >= lower) & (ms1_df["mz"] <= upper)]
            elif level == 2:
                results_df = ms2_df[
                    (ms2_df["precmz"] >= lower) & (ms2_df["precmz"] <= upper)
                ]
            results_df = results_df.copy()
            if not results_df.empty:
                results_df["species"] = species
                results_df["expected_mass"] = expected_mass
                results_df["filename"] = str(relative_filename)
                results_df = results_df.round(
                    {
                        "i": 4,
                        "i_norm": 4,
                        "i_tic_norm": 4,
                        "mz": 4,
                        "rt": 4,
                        "expected_mass": 4,
                    }
                )
                file_path = Path(
                    output_dir,
                    f"{str(row['pubchem_cid'])}_{str(row['metabolite_of_interest'])}_mslevel-{str(level)}.tsv",
                )
                results_df.to_csv(
                    file_path,
                    sep="\t",
                    index=False,
                    mode="a",
                    header=not file_path.exists(),
                )


def process_files_parallel(
    files, metabolites, ion_mass_map, ppm_tolerance, basedir, cpus, level, output_dir
):
    with tqdm(total=len(files), desc="Processing files") as pbar:
        with multiprocessing.Pool(processes=cpus) as pool:
            # Update the progress bar after processing each file
            for _ in pool.imap_unordered(
                partial(
                    process_file,
                    metabolites=metabolites,
                    ion_mass_map=ion_mass_map,
                    ppm_tolerance=ppm_tolerance,
                    basedir=basedir,
                    level=level,
                    output_dir=output_dir,
                ),
                files,
            ):
                pbar.update(1)


def group_final_results(output_dir):
    for file in Path(output_dir).glob("*.tsv"):
        df = pd.read_csv(file, sep="\t")
        df = (
            df.groupby(["filename", "scan"])
            .agg(
                {
                    "i_norm": list,
                    "i_tic_norm": list,
                    "mz": list,
                    "rt": "first",
                    "ms1scan": "first",
                    "precmz": "first",
                    "expected_mass": "first",
                    "charge": "first",
                    "polarity": "first",
                    "species": "first",
                }
            )
            .reset_index()
        )
        # save as feather
        df.to_feather(file.with_suffix(".feather"))


def main(
    metabolites_tsv,
    adducts_tsv,
    basedir,
    subdir,
    ppm_tolerance,
    cpus,
    level,
    output_dir,
    force=False,
):
    if Path(output_dir).exists() and not force:
        print(
            f"Output directory {output_dir} already exists. Use --force to overwrite."
        )
        return
    else:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    metabolites = pd.read_csv(Path(metabolites_tsv), sep="\t")
    metabolites["monoisotopic_mass"] = metabolites["molecular_formula"].apply(
        lambda x: (mass.calculate_mass(formula=x))
    )
    df = pd.read_csv(Path(adducts_tsv), sep="\t")
    ion_mass_map = {
        ion: (mass, df[df["ion"] == ion]["mult"].values[0])
        for ion, mass in zip(df["ion"], df["mass"])
    }
    del df
    inputdir = Path(basedir, subdir)
    files = list(Path(inputdir).glob("**/*.mzXML"))
    process_files_parallel(
        files,
        metabolites,
        ion_mass_map,
        ppm_tolerance,
        basedir,
        cpus,
        level,
        output_dir,
    )
    group_final_results(output_dir)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.level not in [1, 2]:
        raise ValueError("Level must be 1 or 2")
    if args.cpus < 1:
        raise ValueError("Number of cpus must be greater than 0")
    if args.ppm_tolerance < 0:
        raise ValueError("PPM tolerance must be greater than 0")
    main(
        metabolites_tsv=args.metabolites_tsv,
        adducts_tsv=args.adducts_tsv,
        basedir=args.basedir,
        subdir=args.subdir,
        ppm_tolerance=args.ppm_tolerance,
        cpus=args.cpus,
        level=args.level,
        output_dir=args.output_dir,
    )
