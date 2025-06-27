import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from matchms.filtering import normalize_intensities
from matchms.importing import load_from_mzxml


def mirror_plot(
    sample_spectrum,
    ccmslib_spectrum,
    species,
    cos_score,
    n_matches,
    gnps_cluster,
    label_top_n=10,
    max_peaks_per_spectrum=100,
    max_mz=float("inf"),
):
    sample_spectrum = normalize_intensities(sample_spectrum)
    ccmslib_spectrum = normalize_intensities(ccmslib_spectrum)
    sample_df = pd.DataFrame(
        {"mz": sample_spectrum.mz, "intensity": sample_spectrum.intensities}
    )
    ccmslib_df = pd.DataFrame(
        {"mz": ccmslib_spectrum.mz, "intensity": ccmslib_spectrum.intensities}
    )
    sample_df["mz"] = sample_df["mz"].round(4)
    ccmslib_df["mz"] = ccmslib_df["mz"].round(4)
    sample_top_10 = sample_df.nlargest(label_top_n, "intensity")
    ccmslib_top_10 = ccmslib_df.nlargest(label_top_n, "intensity")
    sample_df = sample_df.nlargest(max_peaks_per_spectrum, "intensity")
    ccmslib_df = ccmslib_df.nlargest(max_peaks_per_spectrum, "intensity")
    fig = go.Figure()
    # Add library data with positive peaks as vertical lines (blue)
    for _, row in sample_df.iterrows():
        fig.add_shape(
            type="rect",
            x0=row["mz"] - 0.05,
            x1=row["mz"] + 0.05,  # Horizontal span of the peak
            y0=0,
            y1=row["intensity"],  # The peak stretches from y=0 to its intensity
            line=dict(color="blue", width=2),
            fillcolor="blue",
            opacity=0.5,
        )
    # Add sample data with negative peaks (inverted, red)
    for _, row in ccmslib_df.iterrows():
        fig.add_shape(
            type="rect",
            x0=row["mz"] - 0.05,
            x1=row["mz"] + 0.05,  # Horizontal span of the peak
            y0=0,
            y1=-row[
                "intensity"
            ],  # The peak stretches from y=0 to its negative intensity
            line=dict(color="red", width=2),
            fillcolor="red",
            opacity=0.5,
        )
    # Add label for top peaks in the library (blue)
    fig.add_trace(
        go.Scatter(
            x=sample_top_10["mz"],
            y=sample_top_10["intensity"],
            mode="text",
            text=sample_top_10["mz"].astype(str),
            textposition="top right",
            textfont=dict(size=12, color="blue"),
            name=f"Scan {str(sample_spectrum.metadata.get('scan', None))}",
        )
    )
    # Add label for top peaks in the sample (red)
    fig.add_trace(
        go.Scatter(
            x=ccmslib_top_10["mz"],
            y=-ccmslib_top_10["intensity"],
            mode="text",
            text=ccmslib_top_10["mz"].astype(str),
            textposition="bottom right",
            textfont=dict(size=12, color="red"),
            name=f"GNPS Library ({ccmslib_spectrum.metadata.get('spectrum_id', None)})",
        )
    )
    scan_name = Path(sample_spectrum.metadata.get("file_name", None)).stem
    scan_name = scan_name.removesuffix(".mzXML")
    title = f"Sample {scan_name}; Species: {species}; Score: {cos_score:.2f}; Matches: {n_matches}; GNPS2 Cluster: {gnps_cluster}"
    # Update layout for customization and styling
    fig.update_layout(
        title=title,
        xaxis_title="m/z",
        yaxis_title="Intensity",
        template="plotly_white",  # Light theme
        xaxis_tickangle=90,  # Rotate x-axis labels
        showlegend=True,
    )
    return fig


def plot_gnps_hits(tdf, mgf, GNPS2_RESULTS_DIR, MSV000084900_DIR):
    df = pd.read_csv(
        Path(GNPS2_RESULTS_DIR, "nf_output/clustering/clusterinfo.tsv"), sep="\t"
    )
    tdf = tdf.iloc[[0]]
    first_member = df[df["#ClusterIdx"] == tdf["cluster index"].values[0]].iloc[[3]]
    filename = first_member["#Filename"].values[0]
    filename = Path(filename).stem
    filename = f"{filename}.mzXML"
    filepath = list(MSV000084900_DIR.glob(f"**/*{filename}"))[0]
    spec = None
    for i in load_from_mzxml(str(filepath)):
        if i.metadata["scan_number"] == str(first_member["#Scan"].values[0]):
            spec = i
            break
    temp = spec.metadata
    spec.metadata["file_name"] = filename
    spec.metadata.update({"file_name": filename})
    temp = spec.metadata
    temp.update({"file_name": filename})
    temp.update({"scan": str(first_member["#Scan"].values[0])})
    spec.metadata = temp
    ccmslib = tdf["SpectrumID"].values[0]
    return mirror_plot(
        sample_spectrum=spec,
        ccmslib_spectrum=mgf[ccmslib],
        species=tdf["Adduct"].values[0],
        cos_score=tdf["MQScore"].values[0],
        n_matches=tdf["SharedPeaks"].values[0],
        gnps_cluster=tdf["cluster index"].values[0],
        label_top_n=10,
        max_peaks_per_spectrum=100,
    )
