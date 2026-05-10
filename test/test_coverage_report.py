import pytest

import coverage_report as cov

def test_that_error_raised_if_not_found_file():
    """
    Test that a FileNotFoundError is raised if an input file is not found
    """
    input_file = "sambamba_not_existing.txt"
    with pytest.raises(FileNotFoundError):
        cov.read_sambamba(input_file)