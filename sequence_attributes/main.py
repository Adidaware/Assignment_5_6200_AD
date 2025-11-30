"""
This is the main module of the assignment that imports functions from other modules
This module receives the command line arguments, read the cdds attributes, read the ensemble gene data
"""
import sys
import argparse
import pandas as pd
from sequence_attributes.sequence_formats.fasta_format import get_fasta_lists
from sequence_attributes.utils.seq_attribute_utils import (return_standard_genetic_code,
                                                           get_additional_sequence_attributes)


def main():
    """
    Defining main
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Combine information from multiple files and calculate sequence attributes.")
    parser.add_argument("--infile_ccds_fasta", help="CCDS fna file to open (FASTA format)")
    parser.add_argument("--infile_ccds_attributes", help="CCDS attributes file to open")
    parser.add_argument("--infile_ensembl_gene", help="Ensembl gene information file to open")
    parser.add_argument("--excel_outfile", help="Excel file to write final output")
    args = parser.parse_args()

    # Read CCDS attributes file into pandas dataframe
    ccds_attributes_df = pd.read_csv(args.infile_ccds_attributes, sep="\t")
    ccds_attributes_df.rename(columns={"gene_id": "refseq_gene_id", "#chromosome": "chrom"}, inplace=True)

    # Read Ensembl gene data file into pandas dataframe
    ensembl_gene_df = pd.read_csv(args.infile_ensembl_gene, sep="\t")
    ensembl_gene_df.rename(columns={"ID": "ensembl_gene_id", "Canonical Transcript": "ensembl_canonical_transcript_id",
                                    "Gene": "gene", "Description": "description", "Biotype": "biotype"}, inplace=True)

    # Merge CCDS and Ensembl data frames
    merged_df = pd.merge(ccds_attributes_df, ensembl_gene_df, on="gene", how="left")

    # Get fasta data from the FASTA file
    seq_header, sequence = get_fasta_lists(args.infile_ccds_fasta)

    # Processing the sequence data and storing attributes in a list
    all_data = []
    for i, (header, dna_sequence) in enumerate(zip(seq_header, sequence)):
        all_data.append(get_additional_sequence_attributes(header=header,
                                                           dna_sequence=dna_sequence,
                                                           genetic_code=return_standard_genetic_code(),
                                                           attribute_df=merged_df))
        if i == sys.maxsize:  # Set this to sys.maxsize before final submission
            break

    # Convert list of tuples to pandas dataframe
    result_df = pd.DataFrame(all_data)

    # Sorting the dataframe
    sorted_df = result_df.sort_values(by=["proline_comp", "protein_sequence_len"], ascending=[False, True])

    # Write the sorted dataframe to an Excel file
    excel_outfile = args.excel_outfile.replace(".xlsx", ".tsv")
    sorted_df.to_excel(args.excel_outfile, index=False)

    # Print top 10 and last 10 genes with highest and lowest proline composition
    print("Top 10 genes with highest proline composition:")
    print(sorted_df.head(10)[["ccds", "refseq_gene_id", "biotype", "ensembl_gene_id", "description",
                              "protein_sequence_len", "proline_comp"]])

    print("\nLast 10 genes with lowest proline composition:")
    print(sorted_df.tail(10)[["ccds", "refseq_gene_id", "biotype", "ensembl_gene_id", "description",
                              "protein_sequence_len", "proline_comp"]])

    # Write top and last 10 genes to a Tab-delimited file
    tsv_outfile = excel_outfile.replace(".tsv", "_proline_composition.tsv")
    with open(tsv_outfile, "w", encoding='utf-8') as tsv_file:
        tsv_file.write("Top 10 genes with highest proline composition:\n")
        tsv_file.write(sorted_df.head(10).to_csv(sep="\t", index=False))
        tsv_file.write("\n\nLast 10 genes with lowest proline composition:\n")
        tsv_file.write(sorted_df.tail(10).to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
