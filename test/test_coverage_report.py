import pytest
import pandas as pd

import coverage_report as cov

def test_that_error_raised_if_not_found_file():
    """
    Test that a FileNotFoundError is raised if an input file is not found
    """
    input_file = "sambamba_not_existing.txt"
    with pytest.raises(FileNotFoundError):
        cov.read_sambamba(input_file)


def test_that_missing_gene_detection_is_raised():
    """
    Test that if a gene is not present in the sambamba file, it is printed in the final report as 'MISSING'
    """
    # A sambamba df with ACTA1 data only
    sambamba_df = pd.DataFrame([{
        'chromosome': '1',
        'StartPosition': 1000,
        'FullPosition': '1-1000-2000',
        'EndPosition': 2000,
        'GeneSymbol': 'ACTA1',
        'readCount': 500,
        'meanCoverage': 150.0,
        'percentage30': 100.0,
    }])

    result = cov.create_report(sambamba_df, "panel_genes.txt", 100.0)
    #CFL2 is included in the panel_genes.txt:
    assert result.loc[result['GeneSymbol'] == 'CFL2', 'Status'].values[0] == 'MISSING'