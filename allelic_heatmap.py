#!/usr/bin/env python3

"""
This script can be used to generate a heatmap with the cg/wgMLST alleles distribution across a species tree

By Vitor Borges
@INSA
"""

import sys
import argparse
import datetime as datetime
import pandas as pd
import numpy as np
from Bio import Phylo
import colorsys
from collections import defaultdict
import plotly.graph_objects as go

version = "1.0.0"
last_updated = "2025-09-22"

def lighten_color(hex_color, factor=0.5):
    """
    get fainted colors
    """

    hex_color = hex_color.strip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    h, l, s = colorsys.rgb_to_hls(r/255, g/255, b/255)
    l = min(1, l + (1-l)*factor)
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    
    return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"

def generate_auto_colors(categories):
    """
    get a map of colors
    """
    
    cats = list(categories)
    n = max(1, len(cats))
    color_map_cat = {}
    for i, cat in enumerate(cats):
        h = i / n
        s = 0.6
        v = 0.8
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        color_map_cat[cat] = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
    
    return color_map_cat

def read_tree_order(tree_file):
    """
    get a list with the samples of a tree
    """

    tree = Phylo.read(tree_file, "newick")

    return [term.name for term in tree.get_terminals()]

def assign_colors(allelic_matrix, metadata, meta_col, color_map, mode="frequency", ncat=None, ncat_dominant=None):
    """
    assign colors to the groups
    """
    
    allele_to_color = {}
    exclusives = defaultdict(list)
    shared = defaultdict(list)
    allele_to_categories = defaultdict(set)
    conserved_loci = set()

    all_cats = metadata[meta_col].dropna().unique()
    total_per_cat = metadata[meta_col].value_counts()

    for locus in allelic_matrix.columns:
        locus_data = allelic_matrix[locus]
        locus_data = locus_data.replace("0", np.nan)
        unique_alleles = locus_data.dropna().unique()

        if len(unique_alleles) == 0:
            continue
        elif len(unique_alleles) == 1:
            conserved_loci.add(locus)
            only_allele = unique_alleles[0]
            allele_to_color[(locus, only_allele)] = "#000000"  # conserved = gray
            allele_to_categories[(locus, only_allele)] = set(all_cats)
        else:
            for allele in unique_alleles:
                idx = locus_data[locus_data == allele].index
                categories = metadata.loc[idx, meta_col].dropna()
                counts = categories.value_counts()

                if mode == "frequency":
                    freqs = counts / total_per_cat[counts.index]
                    dominant_cat = freqs.idxmax()
                else:
                    dominant_cat = counts.idxmax()

                allele_to_categories[(locus, allele)] = set(categories.unique())

                # exclusives/shared classification (for tables) → independent of ncat/ncat_dominant
                if len(counts) == 1:
                    base_color = color_map[dominant_cat]
                    allele_to_color[(locus, allele)] = lighten_color(base_color, factor=0.5)
                    exclusives[dominant_cat].append((locus, allele))
                else:
                    base_color = color_map[dominant_cat]
                    allele_to_color[(locus, allele)] = base_color
                    shared[(locus, allele)] = list(categories.unique())

                # Apply "nearly conserved" rules ONLY to color (not to tables)
                if ncat is not None and len(set(categories)) >= ncat:
                    allele_to_color[(locus, allele)] = "#000000"
                if ncat_dominant is not None:
                    # Count in how many categories this allele is dominant
                    dominant_count = 0
                    for cat in all_cats:
                        cat_idx = metadata[metadata[meta_col] == cat].index
                        if len(cat_idx) == 0:
                            continue
                        cat_series = locus_data.loc[cat_idx]
                        if len(cat_series.dropna()) == 0:
                            continue
                        allele_counts = cat_series.value_counts()
                        if allele_counts.idxmax() == allele:
                            dominant_count += 1
                    if dominant_count >= ncat_dominant:
                        allele_to_color[(locus, allele)] = "#000000"

    return allele_to_color, exclusives, shared, allele_to_categories, conserved_loci

def build_heatmap(allelic_matrix, allele_to_color, allele_to_categories, sample_order, color_map, meta_col, output_prefix):
    """
    generate the heatmap
    """

    loci = allelic_matrix.columns.tolist()
    samples = [s for s in sample_order if s in allelic_matrix.index]

    unique_colors = list({v for v in allele_to_color.values()})
    if "#ffffff" not in unique_colors:
        unique_colors.append("#ffffff")
    if "#000000" not in unique_colors:
        unique_colors.append("#000000")
    color_to_id = {color: i for i, color in enumerate(unique_colors)}

    z = []
    hover_text = []
    for sample in samples:
        row = []
        h_row = []
        for locus in loci:
            allele = allelic_matrix.loc[sample, locus]
            if allele == "0" or pd.isna(allele):
                row.append(color_to_id["#ffffff"])  # missing = white
                h_row.append(f"Sample: {sample}<br>{locus}<br>Allele: 0<br>Category: -")
            else:
                color = allele_to_color.get((locus, allele), "#ffffff")
                row.append(color_to_id[color])
                cats = allele_to_categories.get((locus, allele), set())
                if cats:
                    cat_str = ", ".join(sorted(list(cats)))
                else:
                    cat_str = "-"
                h_row.append(f"Sample: {sample}<br>{locus}<br>Allele: {allele}<br>Category: {cat_str}")
        z.append(row)
        hover_text.append(h_row)

    fig = go.Figure(data=go.Heatmap(
        z=z,
        x=[f"Locus {i+1}" for i in range(len(loci))],
        y=samples,
        hoverinfo="text",
        text=hover_text,
        colorscale=[[i/(len(unique_colors)-1), c] for i, c in enumerate(unique_colors)],
        showscale=False
    ))

    fig.update_layout(
        xaxis=dict(showticklabels=False),
        yaxis=dict(showticklabels=False, autorange="reversed"),
        margin=dict(t=50, l=50, r=250),
        title="Allelic Category Heatmap"
    )

    annotations = []
    legend_y = 1
    step = 0.05
    for cat, base_color in color_map.items():
        faint_color = lighten_color(base_color, 0.5)
        annotations.append(dict(x=1.12, y=legend_y, xref='paper', yref='paper',
                                text=f"{cat} (faint=exclusive)", showarrow=False, font=dict(color=faint_color)))
        legend_y -= step
        annotations.append(dict(x=1.12, y=legend_y, xref='paper', yref='paper',
                                text=f"{cat} (dark=shared)", showarrow=False, font=dict(color=base_color)))
        legend_y -= step
    fig.update_layout(annotations=annotations)

    html_out = f"{output_prefix}.html"
    extra_buttons = f"""
<script>
function download(format) {{
    Plotly.downloadImage(document.getElementsByClassName('js-plotly-plot')[0], {{
        format: format,
        filename: '{output_prefix}'
    }})
}}
</script>
<div style="margin-bottom:10px;">
  <button onclick="download('png')">Download PNG</button>
  <button onclick="download('svg')">Download SVG</button>
</div>
"""
    with open(html_out, "w") as f:
        f.write(extra_buttons)
        f.write(fig.to_html(include_plotlyjs='cdn', full_html=False))

    print(f"[OK] Heatmap saved to {html_out}")

def save_summary(exclusives, shared, allele_to_categories, output_prefix):
    """
    generate a summary report
    """

    rows = []
    for cat, alleles in exclusives.items():
        for locus, allele in alleles:
            rows.append([cat, "exclusive", locus, allele])
    for (locus, allele), cats in shared.items():
        cat_str = ",".join(sorted(cats))
        rows.append([cat_str, "shared", locus, allele])
    df = pd.DataFrame(rows, columns=["Category", "Type", "Locus", "Allele"])
    out1 = f"{output_prefix}_allele_summary.tsv"
    df.to_csv(out1, sep="\t", index=False)
    print(f"[OK] Allele summary saved to {out1}")

    # Summary counts (independent of ncat/ncat_dominant)
    summary_counts = (
        df.groupby(["Category", "Type"])
        .agg(
            n_loci=("Locus", "nunique"),
            n_alleles=("Allele", "count")  # total occurrences across loci
        )
        .reset_index()
    )
    out2 = f"{output_prefix}_allele_summary_counts.tsv"
    summary_counts.to_csv(out2, sep="\t", index=False)
    print(f"[OK] Allele summary counts saved to {out2}")

def save_auto_colors_tsv(color_map, output_prefix, meta_col):
    """
    generate tsv with the color schema selected
    """

    rows = []
    for cat, color in color_map.items():
        rows.append([meta_col, cat, color])
    df = pd.DataFrame(rows, columns=["column", "category", "color"])
    out = f"{output_prefix}_colors.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[OK] Auto-generated colors saved to {out}")

def build_detailed_matrix(allelic_matrix, sample_order, allele_to_categories, conserved_loci, output_prefix):
    """
    generate a detailed matrix with the allele presence/absence
    """

    loci = allelic_matrix.columns.tolist()
    samples = [s for s in sample_order if s in allelic_matrix.index]

    data = []
    for sample in samples:
        row = []
        for locus in loci:
            allele = allelic_matrix.loc[sample, locus]
            if allele == "0" or pd.isna(allele):
                row.append("0 (missing; -)")
            else:
                cats = allele_to_categories.get((locus, allele), set())
                if locus in conserved_loci:
                    annot = f"{allele} (conserved; all)"
                else:
                    if len(cats) <= 1:
                        cat_str = "-" if len(cats) == 0 else ", ".join(sorted(list(cats)))
                        annot = f"{allele} (exclusive; {cat_str})"
                    else:
                        cat_str = ", ".join(sorted(list(cats)))
                        annot = f"{allele} (shared; {cat_str})"
                row.append(annot)
        data.append(row)

    df = pd.DataFrame(data, index=samples, columns=loci)
    out = f"{output_prefix}_allelic_detailed_matrix.tsv"
    df.to_csv(out, sep="\t")
    print(f"[OK] Detailed matrix saved to {out}")

def build_dominance_table(allelic_matrix, metadata, meta_col, mode, output_prefix):
    """
    generate a matrix with information about the dominant alleles
    """

    all_cats = list(metadata[meta_col].dropna().unique())
    total_per_cat = metadata[meta_col].value_counts()
    rows = []

    for locus in allelic_matrix.columns:
        locus_series = allelic_matrix[locus]
        present = locus_series.replace("0", np.nan).dropna()
        unique_alleles = present.unique()
        for allele in unique_alleles:
            idx = present[present == allele].index
            cats = metadata.loc[idx, meta_col].dropna()
            counts = cats.value_counts()

            count_map = {cat: int(counts.get(cat, 0)) for cat in all_cats}
            freq_map = {cat: (count_map[cat] / total_per_cat[cat]) if cat in total_per_cat else 0.0 for cat in all_cats}

            if mode == "frequency":
                dominant_cat = max(freq_map, key=freq_map.get)
            else:
                dominant_cat = max(count_map, key=count_map.get)

            type_val = "exclusive" if sum(1 for v in count_map.values() if v > 0) == 1 else "shared"

            counts_str = "; ".join([f"{cat}:{count_map[cat]}" for cat in sorted(all_cats)])
            freqs_str = "; ".join([f"{cat}:{freq_map[cat]:.4f}" for cat in sorted(all_cats)])

            rows.append([locus, allele, dominant_cat, type_val, counts_str, freqs_str])

    df = pd.DataFrame(rows, columns=["Locus", "Allele", "Dominant_Category", "Type", "Counts_per_category", "Frequencies_per_category"])
    out = f"{output_prefix}_allele_dominance.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[OK] Allele dominance table saved to {out}")

def main():
    parser = argparse.ArgumentParser(description="Fast allelic heatmap ordered by tree")
    parser.add_argument("-a", "--allelic_matrix", required=True, help="TSV allelic matrix")
    parser.add_argument("-m", "--metadata", required=True, help="Metadata TSV")
    parser.add_argument("-c", "--colors", required=False, help="(Optional) Colors TSV (column, category, color). If omitted, colors are auto-assigned and saved to <prefix>_colors.tsv")
    parser.add_argument("-hcol", "--header", required=True, help="Metadata column to use for categories")
    parser.add_argument("-t", "--tree", required=True, help="Newick tree file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix")
    parser.add_argument("-mode", "--mode", choices=["frequency", "count"], default="frequency", help="Dominant category mode: 'frequency' (default) or 'count'")
    parser.add_argument("-ncat", "--ncat", type=int, default=None, help="Optional: if an allele is present in at least N categories, color it gray (nearly conserved)")
    parser.add_argument("-ncat_dominant", "--ncat_dominant", type=int, default=None, help="Optional: if an allele is dominant in at least N categories, color it gray (nearly conserved)")
    
    args = parser.parse_args()

    print("\n******************** running allelic_heatmap.py ********************\n")
    print("version " + str(version) + " last updated on " + str(last_updated) + "\n")
    print(" ".join(sys.argv))
	
    start = datetime.datetime.now()
    print("start: " + str(start))

    allelic_matrix = pd.read_csv(args.allelic_matrix, sep="\t", index_col=0, dtype=str)
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=0, dtype=str)
    sample_order = read_tree_order(args.tree)

    if args.colors:
        colors_df = pd.read_csv(args.colors, sep="\t", dtype=str)
        csub = colors_df[colors_df["column"] == args.header]
        if csub.empty:
            raise ValueError(f"No colors found in {args.colors} for column '{args.header}'.")
        color_map = {row["category"]: row["color"] for _, row in csub.iterrows()}
    else:
        cats = metadata[args.header].dropna().unique()
        auto = generate_auto_colors(cats)
        color_map = auto
        save_auto_colors_tsv(color_map, args.output_prefix, args.header)

    allele_to_color, exclusives, shared, allele_to_categories, conserved_loci = assign_colors(
        allelic_matrix, metadata, args.header, color_map, args.mode, args.ncat, args.ncat_dominant
    )

    build_heatmap(allelic_matrix, allele_to_color, allele_to_categories, sample_order, color_map, args.header, args.output_prefix)
    save_summary(exclusives, shared, allele_to_categories, args.output_prefix)
    build_detailed_matrix(allelic_matrix, sample_order, allele_to_categories, conserved_loci, args.output_prefix)
    build_dominance_table(allelic_matrix, metadata, args.header, args.mode, args.output_prefix)

    end = datetime.datetime.now()
	
    print("end: " + str(end))
    print("\nDone!")
    print("************************************************************")

if __name__ == "__main__":
    main()