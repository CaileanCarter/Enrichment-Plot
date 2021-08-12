# !/usr/bin/python3

# Plot a Cluster Map of gene sets from Roary and BioCyc's SmartTable of Enriched Pathways.
import argparse
import glob
import re
from os import path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser()

parser.add_argument("RoaryFile", metavar="file path", type=str, help="File path for Roary gene presence/absence data")

parser.add_argument("-s", "--SmartTable", type=str, metavar="file path", help="File or folder path for SmartTables generated from Roary data for organism databases of your choice")

parser.add_argument("-a", "--accessories", action="store_true", default=False, help="Show accessory pathways only.")

parser.add_argument("-c", "--core", action="store_true", default=False, help="Print core pathways.")

parser.add_argument("--score", type=float, metavar="0 <= score <= 1", default=1, help="Assign cutoff score for whether a pathway or core of accessory")

parser.add_argument("--output", type=str, metavar=r"\file\path\filename.xslx", help="Output pathway frequency to an Excel file")

args = parser.parse_args()


relpath = path.dirname(path.abspath(__file__))

def open_roary(RoaryData):
    print("Reading Roary gene presence/absence file...")
    GPA = pd.read_table(RoaryData, index_col=0)
    GPAfilter = generate_clean_gene_set(GPA)

    return GPAfilter, GPA.columns


def generate_clean_gene_set(GPA):
    GPA_filter = {}
    for row in GPA.itertuples():
        search = re.split(r"_\d+|\d+$", row[0])
        gene = search[0]
        if gene in GPA_filter.keys():
            GPA_filter[gene] = np.maximum(GPA_filter[gene], row[1:])
        else:
            GPA_filter[gene] = row[1:]

    return GPA_filter


def output_clean_gene_set(GPA_filter):
    with open(path.join(f"{relpath}" + "roarygenes.txt"), "w") as output:
        for geneID in GPA_filter.keys():
            output.write(geneID + "\n")
    print(f"Common genes outputted to: \n{relpath}\\roarygenes.txt")


def PlotAll(df):
    sns.clustermap(df, cmap="mako")                                           # Plot results
    print("Done.")
    plt.show()


def PlotAccessories(df, cutoff):
    accessory = df[df.sum(axis=1) < cutoff]
    sns.clustermap(accessory, cmap="mako")                                    # Plot results
    print("Done")
    plt.show()


def PrintCore(df, cutoff):
    core = df[df.sum(axis=1) >= cutoff]
    print("Core pathways:")
    for pathway in core.index:
        print(pathway)


def Single_enrichment_file(filepath):
    enrichment = pd.read_table(filepath, index_col=0)                         # Load BioCyc's SmartTable Enrichment table
    enrichment = enrichment["Matches"].str.split(" // ", expand=True)
                                                   
    enrichdict = enrichment.to_dict("index")                                  # Convert enrichment table to dictionary (easier to work with)
    
    filtered_dict = {Pathway : [gene for gene in genes.values() if gene] 
        for Pathway, genes in enrichdict.items()}                             # Remove all empty data (Pandas has no way to remove empty sets during conversion)

    print("File OKAY.")
    return filtered_dict


def Multi_enrichment_files(folder):
    sharedpath = rf"{folder}\**\\"
    allenrichments = glob.glob(path.join(sharedpath, "Enrich*.txt"))                                    # make a list of files
    enrichments = pd.concat((pd.read_csv(file, sep="\t", index_col=0) for file in allenrichments))      # combine all enrichment data
    enrichments = enrichments["Matches"].str.split(" // ", expand=True)                                 # split the matches into own cells
    
    enrichdict = {}
    for row in enrichments.itertuples():                                                                # convert into a dictionary so duplicate pathways have their genes merged
        pathway = row[0]
        if pathway not in enrichdict.keys():
            enrichdict[pathway] = [gene for gene in row[1:] if gene] 
        else:
            addlist = [gene for gene in row[1:] if gene]
            enrichdict[pathway] += addlist

    enrichdict_toset = {pathway : set(genes) for pathway, genes in enrichdict.items()}                  # convert to a set
    print("Enrichment files have merged successfully.")
    return enrichdict_toset


def Main(RoaryData, EnrichmentData, accessories=False, core=False, score=1, output=None):

    if path.isdir(EnrichmentData):
        print("Folder path given, attempting to merge enrichment files...")
        enrichmentdict = Multi_enrichment_files(EnrichmentData)
    elif path.isfile(EnrichmentData):
        print("File path given, reading file...")
        enrichmentdict = Single_enrichment_file(EnrichmentData)
    else:
        raise IOError(f"No file or folder path given. {EnrichmentData} is not a valid path.")
    
    GPA_filter, columns = open_roary(RoaryData)
    condensedGPA = pd.DataFrame.from_dict(GPA_filter, orient="index", columns=columns)
    print("Done.")

    print("Allocating pathways to isolates...")
    PathwayDict = {}                                                                  # Keep count of how many times a gene occurs for each pathway for each isolate ID
    for Pathway, genes in enrichmentdict.items():
        pathwaycount = {isolateID : 0 for isolateID in condensedGPA.columns}
        for gene in genes:
            try:
                gene_presence_absence = condensedGPA.loc[gene]                        # Pull gene presence/absence from Roary data
            except KeyError:
                continue

            for isolateID, occur in gene_presence_absence.to_dict().items():
                pathwaycount[isolateID] += occur
        PathwayDict[Pathway] = pathwaycount    

    PathwayTable = pd.DataFrame.from_dict(PathwayDict, orient="index")                # Transform to DataFrame for plotting
    PathwayTable = PathwayTable.apply(lambda row: row / row.max(), axis=1)            # normalise occurance into frequency

    if 0 <= score <= 1:
        cutoff = len(PathwayTable.columns) * score
    else:
        raise ValueError(f"Score must be 0 <= score <= 1, not {score}")

    if output:
        print(f"Creating Excel file at location: {output}")
        PathwayTable.to_excel(output)
    if core:
        print("Printing core pathways...")
        PrintCore(PathwayTable, cutoff)
    if accessories:
        print("Plotting accessory pathways...")
        PlotAccessories(PathwayTable, cutoff)
    else:
        print("Plotting all pathways..")
        PlotAll(PathwayTable)
    

if __name__ == "__main__":

    if not args.SmartTable:
        GPA, _ = open_roary(args.RoaryFile)
        output_clean_gene_set(GPA)

    else:
        Main(args.RoaryFile,
            args.SmartTable, 
            accessories=args.accessories, 
            core=args.core, 
            score=args.score, 
            output=args.output
            )


