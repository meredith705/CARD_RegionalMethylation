import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import gzip
import sys

"""
python plotMethylViolin.py.py regionsToPlot.bed large_sample_data.bed violin_plot.png

"""

def get_header(path, gzipped=False):
    opener = gzip.open if gzipped else open
    with opener(path, 'rt') as f:
        return f.readline().strip().split('\t')

def bedtool_without_header(bed_gz_path):
    # Return a BedTool object that skips the header line from a bgzipped file
    def line_gen():
        with gzip.open(bed_gz_path, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue  # skip header
                yield line
    return pybedtools.BedTool(line_gen())

def main(small_bed_path, large_bed_path, output_plot):

    # store headers
    small_header = get_header(small_bed_path)
    large_header = get_header(large_bed_path, gzipped=True)

    # Load BEDs using pybedtools
    small_bed = pybedtools.BedTool(small_bed_path)
    # large_bed = pybedtools.BedTool(large_bed_path)
    large_bed = bedtool_without_header(large_bed_path)

    # Intersect group bed with regions to be plotted
    intersected = large_bed.intersect(small_bed, wa=True, wb=True)

    if len(intersected) == 0:
        print("No intersections found!")
        sys.exit(1)

    # Column names for intersected file:
    # [group fields] + [plotting bed fields]
    column_names = large_header + [f"{h}_b" for h in small_header]

    # Convert to DataFrame with known column names
    df = intersected.to_dataframe(names=column_names, header=None)

    print(f'intersected data {df.shape[0]} regions to plot:\n',df.head())

    # Extract region labels from small bed (which now has suffix _b)
    df['region'] = df[[f'{small_header[0]}_b', f'{small_header[1]}_b', f'{small_header[2]}_b']].astype(str).agg(':'.join, axis=1)

    # Exdlude columns from the input bed - which should be present at the beginning and the end of the intersecton
    columns_to_exclude_from_plot = len(small_header)
    sample_cols = large_header[columns_to_exclude_from_plot:-columns_to_exclude_from_plot]  
    # Melt into long format for seaborn plotting
    melted = df.melt(id_vars=['region'], value_vars=sample_cols,
                     var_name='sample', value_name='AverageMethylation')
    melted['AverageMethylation'] = pd.to_numeric(melted['AverageMethylation'], errors='coerce')

    print('melted plot data:\n',melted.head())


    # Plot
    # Adjust width to compensate for the number of regions plotted
    figwidth = 4 * df.shape[0]
    plt.figure(figsize=(figwidth, 8))

    ax = sns.violinplot(
        x='region',
        y='AverageMethylation',
        data=melted,
        inner='box',
        # density_norm='width',
        alpha=0.75,              
        linewidth=1,
        palette='husl',
        hue='region',
        legend=False
    )
    sns.swarmplot(
        x='region',
        y='AverageMethylation',
        data=melted,
        color='black',
        size=2,
        ax=ax
    )

    plt.xticks(rotation=45, ha='right')

    plt.ylim(-5, 105)
    
    plt.title("Average Methylation per Region")

    plt.tight_layout()

    plt.savefig(output_plot)
    print(f"Saved violin plot to: {output_plot}")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python bed_violin_plot.py <input_regions.bed> <large_data.bed> <output_plot.png>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
