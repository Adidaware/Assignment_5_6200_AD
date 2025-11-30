"""
Unit tests for io_utils.py module
"""
import pytest
from sequence_attributes.utils.io_utils import FileHandler


def test_file_handler_success():
    """
    Test that FileHandler successfully opens a file and closes it.
    """
    with FileHandler('test_file.txt', mode='w') as _f:
        _f.write("Test content")
    assert _f.closed


def test_file_handler_open_failure():
    """
    Test that FileHandler raises an OSError when it fails to open a file.
    """
    with pytest.raises(OSError):
        with FileHandler('non_existent_file.txt', mode='r'):
            pass


def test_file_handler_invalid_mode():
    """
    Test that FileHandler raises a ValueError when an invalid mode is provided.
    """
    with pytest.raises(ValueError):
        with FileHandler('test_file.txt', mode='invalid_mode'):
            pass


def test_file_handler_invalid_type():
    """
    Test that FileHandler raises a TypeError when an invalid type is provided for mode.
    """
    with pytest.raises(TypeError):
        with FileHandler('test_file.txt', mode=123):
            pass


def test_file_handler_close_on_exception():
    """
    Test that FileHandler closes the file even if an exception is raised within the with block.
    """
    try:
        with FileHandler('test_file.txt', mode='w') as _f:
            raise IOError("Test exception")
    except IOError:
        pass
    assert _f.closed
