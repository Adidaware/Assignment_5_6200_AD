"""
This module contains all the functions which are imported in main.py
"""
from typing import List, Union
from collections import namedtuple
import pandas as pd


def gc_content(sequence: str, round_to: int = 2, percentage: bool = True) -> float:
    """ Calculate the GC content of a DNA sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    gc_percentage = (gc_count / total_bases) * 100 if percentage else gc_count / total_bases
    return round(gc_percentage, round_to)


def lookup_by_ccds(ccds_id: str, chrom: str, df: pd.DataFrame) -> pd.DataFrame:
    """Lookup CCDS entry by CCDS ID and chromosome."""
    return df[(df['ccds_id'] == ccds_id) & (df['chrom'] == chrom)]


def get_tm_from_dna_sequence(_sequence: str, _round_to: int = 2) -> float:
    """Calculate the melting temperature (Tm) for a DNA sequence."""
    # Implementation of Tm calculation here
    return 0.0  # Placeholder


def get_sequence_composition(sequence: Union[str, list, tuple], sequence_type: str = 'dna') -> dict:
    """
    Compute the composition of a DNA or protein sequence.

    @param sequence: The sequence to compute the composition of.
    @param sequence_type: The type of sequence. Can be either 'dna' or 'protein'.
    @return: A dictionary containing the composition of the sequence.
    """
    if sequence_type == 'dna':
        composition = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for base in sequence:
            composition[base] += 1
    elif sequence_type == 'protein':
        composition = {amino_acid: sequence.count(amino_acid) for amino_acid in set(sequence)}
    else:
        raise ValueError("Invalid sequence type. Must be either 'dna' or 'protein'.")
    return composition


def calculate_amino_acid_content(amino_acid: str, sequence: str, round_to: int = 2, percentage: bool = True) -> float:
    """Calculate the content of a specific amino acid in a protein sequence."""
    aa_count = sequence.count(amino_acid)
    total_aa = len(sequence)
    if total_aa == 0:
        return 0 if percentage else 0.0
    aa_percentage = (aa_count / total_aa) * 100 if percentage else aa_count / total_aa
    return round(aa_percentage, round_to)


def extract_kmers(sequence: str, k: int = 3) -> List[str]:
    """Extract k-mers from a DNA sequence."""
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmers.append(sequence[i:i + k])
    return kmers


def return_standard_genetic_code() -> dict:
    """Return the standard genetic code."""
    return {
        "UUU": "F", "UUC": "F", "UCU": "S", "UCC": "S", "UAU": "Y", "UAC": "Y",
        "UGU": "C", "UGC": "C", "UUA": "L", "UCA": "S", "UAA": None, "UGA": None,
        "UUG": "L", "UCG": "S", "UAG": None, "UGG": "W", "CUU": "L", "CUC": "L",
        "CCU": "P", "CCC": "P", "CAU": "H", "CAC": "H", "CGU": "R", "CGC": "R",
        "CUA": "L", "CUG": "L", "CCA": "P", "CCG": "P", "CAA": "Q", "CAG": "Q",
        "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "ACU": "T", "ACC": "T",
        "AAU": "N", "AAC": "N", "AGU": "S", "AGC": "S", "AUA": "I", "ACA": "T",
        "AAA": "K", "AGA": "R", "AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
        "GUU": "V", "GUC": "V", "GCU": "A", "GCC": "A", "GAU": "D", "GAC": "D",
        "GGU": "G", "GGC": "G", "GUA": "V", "GUG": "V", "GCA": "A", "GCG": "A",
        "GAA": "E", "GAG": "E", "GGA": "G", "GGG": "G"
    }


def protein_translation(sequence: str, genetic_code: dict) -> str:
    """
    Translate a DNA sequence into a protein sequence using a given genetic code.

    Parameters:
    - sequence: A string representing the DNA sequence.
    - genetic_code: A dictionary mapping 3-letter RNA codons to 1-letter amino acids.

    Returns:
    - A string representing the translated protein sequence.
    """
    # Initialize the protein sequence
    protein_sequence = ''

    # Iterate over the DNA sequence in steps of 3
    for i in range(0, len(sequence), 3):
        # Extract the codon from the sequence
        codon = sequence[i:i + 3]

        # Convert the codon to RNA by replacing 'T' with 'U'
        rna_codon = codon.replace('T', 'U')

        # Translate the RNA codon to an amino acid using the genetic code
        amino_acid = genetic_code.get(rna_codon, None)

        # If the amino acid is None (stop codon), stop the translation
        if amino_acid is None:
            break

        # Append the amino acid to the protein sequence
        protein_sequence += amino_acid

    return protein_sequence


def get_additional_sequence_attributes(header: str = None,
                                       dna_sequence: str = None,
                                       genetic_code: dict = None,
                                       attribute_df: pd.DataFrame = None) -> namedtuple:
    """
    Get additional attributes for the Header and DNA sequence passed in
    @param header: FASTA header
    @param dna_sequence: FASTA sequence
    @param genetic_code: dictionary where key = codon and key = amino acid
    @param attribute_df: The attribute data frame (gene level information)
    @return: namedtuple
    """
    ccds_id, _, chrom = header.split("|")
    chrom = chrom.replace('chr', '')
    ccds_entry_df = lookup_by_ccds(ccds_id=ccds_id, chrom=chrom, df=attribute_df)

    # Define attributes for the namedtuple
    ccds_attributes_to_get = ['nc_accession', 'gene', 'refseq_gene_id', 'biotype', 'ensembl_gene_id',
                              'ensembl_canonical_transcript_id', 'description']
    tuple_fields = ['ccds', 'chrom', 'header', 'dna_sequence', 'dna_sequence_len', 'protein_sequence',
                    'protein_sequence_len', 'nucleotide_comp', 'amino_acid_comp', 'kmers_comp',
                    'proline_comp', 'tm', 'gc'] + ccds_attributes_to_get

    # Create the named tuple type
    FastaTuple = namedtuple("FastaTuple", tuple_fields)

    # 1. Extract 3-mers from the DNA sequence
    kmers = extract_kmers(dna_sequence)

    # 2. Translate the DNA sequence to a protein sequence
    protein_sequence = protein_translation(dna_sequence, genetic_code)

    # 3. Get the TM of the DNA sequence
    tm_val = get_tm_from_dna_sequence(dna_sequence)

    # 4. Compute compositions for nt, amino acids, and kmers
    nucleotide_comp = get_sequence_composition(dna_sequence, sequence_type='dna')
    amino_acid_comp = get_sequence_composition(protein_sequence, sequence_type='protein') if protein_sequence else None
    kmers_comp = {kmer: kmers.count(kmer) for kmer in set(kmers)}

    # 5. Get GC content
    gc_val = gc_content(dna_sequence)

    # 6. Get CCDS attributes from the DF
    ccds_attributes = [ccds_entry_df[val].iloc[0] for val in ccds_attributes_to_get]

    # 7. Find proline composition
    proline_comp = calculate_amino_acid_content('P', protein_sequence) if protein_sequence else None

    # Instantiate and return the named tuple
    return FastaTuple(ccds_id, chrom, header, dna_sequence, len(dna_sequence), protein_sequence,
                      len(protein_sequence), nucleotide_comp, amino_acid_comp, kmers_comp,
                      proline_comp, tm_val, gc_val, *ccds_attributes)


def _example_sequence():
    """Example DNA sequence for testing."""
    return "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"


def _example_genetic_code():
    """Example genetic code for testing."""
    return return_standard_genetic_code()
