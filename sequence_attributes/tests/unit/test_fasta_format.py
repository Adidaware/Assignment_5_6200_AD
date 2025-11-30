"""
Testing code for fasta_format.py
"""
import os
import pytest
from sequence_attributes.sequence_formats.fasta_format import get_fasta_lists, _verify_lists


@pytest.fixture
def _sample_fasta(tmpdir):
    """
    Defining function for sample FASTA
    """
    fasta_content = ">seq1\nATCGATCG\n>seq2\nGATCGATC"
    file_path = os.path.join(tmpdir, "sample.fasta")
    with open(file_path, "w", encoding="utf-8") as _f:
        _f.write(fasta_content)
    return file_path


def test_get_fasta_lists(_sample_fasta):
    """
    Testing function for get_fasta_lists
    """
    headers, sequences = get_fasta_lists(_sample_fasta)
    assert len(headers) == 2
    assert len(sequences) == 2
    assert headers == ["seq1", "seq2"]
    assert sequences == ["ATCGATCG", "GATCGATC"]


def test_get_fasta_lists_empty_file(tmpdir):
    """
    Testing function for get_fasta_lists
    """
    empty_file_path = os.path.join(tmpdir, "empty.fasta")
    with open(empty_file_path, "w", encoding="utf-8"):
        pass  # Create an empty file
    headers, sequences = get_fasta_lists(empty_file_path)
    assert len(headers) == 0
    assert len(sequences) == 0


def test_get_fasta_lists_invalid_file():
    """
    Testing function for get_fasta_lists if the file does not exist
    """
    with pytest.raises(Exception):
        get_fasta_lists("non_existent_file.fasta")


def test_verify_lists():
    """Test the _verify_lists function."""
    # Test with two lists of the same length
    assert _verify_lists(['header1', 'header2'], ['seq1', 'seq2']) == True, "Same length lists should return True"

    # Test with two lists of different lengths
    assert _verify_lists(['header1', 'header2'], ['seq1']) == False, "Different length lists should return False"

    # Test with one empty list and one non-empty list
    assert _verify_lists(['header1'], []) == False, "One empty list should return False"

    # Test with both lists being empty
    assert _verify_lists([], []) == True, "Both empty lists should return True"

    # Test with one list being None and the other being a list
    assert _verify_lists(None, ['seq1']) == False, "One None list should return False"

    # Test with both lists being None
    assert _verify_lists(None, None) == True, "Both None lists should return True"
