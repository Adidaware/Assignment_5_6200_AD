"""
This module contains the functions generate a list of headers and sequences
"""
from typing import Tuple
from sequence_attributes.utils.io_utils import FileHandler


def get_fasta_lists(infile: str = None) -> Tuple[list, list]:
    """Get sequences and headers from a FASTA file."""
    header_list = []
    seq_list = []
    with FileHandler(infile, mode='r', encoding='utf-8') as fh_in:
        header = ''
        sequence = ''
        for line in fh_in:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence:
                    header_list.append(header)
                    seq_list.append(sequence)
                    sequence = ''
                header = line[1:]
            else:
                sequence += line
        # Add the last sequence
        if header and sequence:
            header_list.append(header)
            seq_list.append(sequence)
    # Verify lists
    if not _verify_lists(header_list, seq_list):
        raise ValueError("Header and sequence lists have different lengths.")
    return header_list, seq_list


def _verify_lists(header_list: list = None, seq_list: list = None) -> bool:
    """Verify if the sizes of header and sequence lists are the same."""
    # Check if either list is None and handle accordingly
    if header_list is None or seq_list is None:
        if header_list is None and seq_list is None:
            # If both are None, it's a valid case (assuming you want to treat this as equal)
            return True
        else:
            # If only one is None, they are not of the same length
            print("One of the lists is None.")
            return False

    # Now that we've handled None cases, we can safely compare the lengths
    if len(header_list) != len(seq_list):
        print("Header and sequence lists have different lengths.")
        return False
    return True
