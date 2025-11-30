"""
Unit tests for seq_attributes_utils.py module
"""
import pytest
import pandas as pd
from sequence_attributes.utils.seq_attribute_utils import (
    gc_content, lookup_by_ccds, get_tm_from_dna_sequence, get_sequence_composition,
    calculate_amino_acid_content, extract_kmers, return_standard_genetic_code, protein_translation,
    get_additional_sequence_attributes, _example_genetic_code)


def test_gc_content():
    """Test gc_content function."""
    sequence = "ATGCGATCGATCG"
    assert gc_content(sequence) == pytest.approx(53.85)


def test_lookup_by_ccds():
    """Test lookup_by_ccds function."""
    df = pd.DataFrame({'ccds_id': ['CCDS1', 'CCDS2'], 'chrom': ['1', '2'], 'value': [10, 20]})
    result = lookup_by_ccds('CCDS1', '1', df)
    assert len(result) == 1
    assert result['value'].iloc[0] == 10


def test_get_tm_from_dna_sequence():
    """Test get_tm_from_dna_sequence function."""
    sequence = "ATGCGATCGATCG"
    assert get_tm_from_dna_sequence(sequence) == 0.0


def test_get_sequence_composition():
    """Test get_sequence_composition function."""
    sequence = "ATGCGATCGATCG"
    composition = get_sequence_composition(sequence)
    assert composition == {'A': 3, 'T': 3, 'G': 4, 'C': 3}


def test_calculate_amino_acid_content():
    """Test calculate_amino_acid_content function."""
    assert calculate_amino_acid_content('P', 'APPPAP') == 66.67


def test_extract_kmers():
    """Test extract_kmers function."""
    assert extract_kmers('ATGCGA', 3) == ['ATG', 'TGC', 'GCG', 'CGA']


def test_return_standard_genetic_code():
    """Test return_standard_genetic_code function."""
    genetic_code = return_standard_genetic_code()
    assert len(genetic_code) == 64
    assert genetic_code["UUU"] == "F"


def test_protein_translation():
    """Test protein_translation function."""
    sequence = "ATG"
    genetic_code = return_standard_genetic_code()
    assert protein_translation(sequence, genetic_code) == "M"


def test_get_additional_sequence_attributes_no_protein():
    """
    Test function for defined function 'get_additional_sequence_attributes'
    """
    header = "CCDS1|1|1"
    dna_sequence = "ATGCGCAT"
    genetic_code = _example_genetic_code()
    attribute_df = pd.DataFrame(
        {"ccds_id": ["CCDS1"], "chrom": ["1"], "nc_accession": ["NM_1"], "gene": ["GENE1"], "refseq_gene_id": ["1"],
         "biotype": ["protein_coding"], "ensembl_gene_id": ["ENSG1"],
         "ensembl_canonical_transcript_id": ["ENST1"], "description": ["Gene 1"]})
    result = get_additional_sequence_attributes(header, dna_sequence, genetic_code, attribute_df)
    assert result.protein_sequence
