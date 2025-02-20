# Assignment_5_6200_AD

File name: main.py
Author's Name: Aditya Daware
File description:

Modules used:
sys
argparse
pandas
sequence_attributes.sequence_formats.fasta_format
sequence_attributes.utils.seq_attribute_utils

Functions used in this program 
: def 
: with open 
: / 
: print 
: if, else 
: for
: argv.parser

Code function: This is the main module of the assignment that imports functions from other modules
This module receives the command line arguments, read the cdds attributes, read the ensemble gene data

File name: io_utils.py
File description: This code is provided in the assignment.

Modules used:
sys

Code function: This module utilizes the context management protocol (implementing __enter__ and __exit__ methods) and presents a robust approach for managing file operations within a large project


File name: seq_attribute_utils.py
Author's Name: Aditya Daware
File description:

Modules imported: 
typing
collections
pandas

Functions used in this program 
: def 
: len() 
: / 
: if, else 
: try, except
: for
: range
: return

Code function: This code is designed to contain all the functions required to read and process cdds data, ensemble gene data and generate the header list, sequence list. All the functions are imported in the main.py module


File name: fasta_format.py
Author's Name: Aditya Daware
File description:

Modules imported: 
typing
sequence_atrributes.utils.io_utils

Functions used in this program 
: def 
: / 
: if, else 
: try, except 
: return
: len()
: ValuError

Code function: This module contains the functions generate a list of headers and sequences


File name: test_fasta_format.py
Author's Name: Aditya Daware

File description:

Modules imported: 
pytest
sequence_attributes.sequence_formats.fasta_format

Functions used in this program 
: def 
: len() 
: .join
: print 
: if, else 
: with open
: assert
: pass

Code function: This python program is designed to test the functions defined in the fasta_format.py module


File name: test_io_utils.py
Author's Name: Aditya Daware
File description:

Modules imported: 
pytest
sequence_attributes.utils.io_utils

Functions used in this program 
: def 
: assert 
: with
: ValueError
: pass
: try, except

Code function: This python program is designed to test the functions in io_utils.py


File name: test_seq_attribute_utils.py
Author's Name: Aditya Daware
File description:

Modules imported: 
pytest
pandas
sequence_attributes.utils.seq_attribute_utils.py

Functions used in this program 
: def 
: assert
: len()
: with
: ValueError
: pass
: try, except

Code function: This python program is designed to test the functions defined in seq_attribute_utils.py