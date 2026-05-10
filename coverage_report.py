#Rationale: process a sambamba output file and create a coverage report highlighting any genes that do not meet the 100% at 30× threshold.

# system
import argparse
import os

#library from 3rd parties:
import pandas as pd

def parse_args() -> argparse.Namespace:
    """
    Parse the inputs given on the command line
    
    Returns
    ----------
    args : Namespace
        Namespace object of passed command line argument inputs
    """    
    parser = argparse.ArgumentParser(
        description="Calculate gene with less than 100% 30x coverage from a sambamba_output.txt file"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="The sambamba_output.txt file to be processed"
    )

    parser.add_argument(
        "-g",
        "--genelist",
        type=str,
        required=True,
        help="The gene list file"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=True,
        help="The coverage threshold"
    )     
    
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="The report.txt file with the coverage report"
    )
    
    args = parser.parse_args()
    return args

def read_sambamba(input_file: str) -> pd.DataFrame:
    """
    Read a txt file to a Pandas dataframe, check for any anomalous rows, handles the gene name column
    #The input may contain rows with missing gene symbols or non-numeric values in the percentage30 column.

    Parameters
    ----------
    input_file : str
        Name of input file

    Returns
    ----------
    pd.DataFrame
        txt file contents as Pandas dataframe with 'GeneSymbol' and 'Accession' separated; filtered for essential rows
    
    Raises
    ----------
    SystemExit
        If file not found
    """

    print("---Parse the sambamba file---")

    try:
        dataframe = pd.read_csv(input_file, sep=r'\s+', header=0)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"File not found -> {input_file}.") from exc
    
    
    # remove all special characters using regex
    dataframe.columns = (dataframe.columns.str.replace(r"[^\w]", "", regex = True))

    #divide the GeneSymbol;Accession column
    dataframe[['GeneSymbol', 'Accession']] = dataframe["GeneSymbolAccession"].str.split(';', n=1, expand=True)

    #keep columns:
    dataframe = dataframe[['chromosome','StartPosition', 'FullPosition','EndPosition', 'GeneSymbol', 'readCount', 'meanCoverage', 'percentage30']]
    
    #Drop duplicates in case there are:
    dataframe = dataframe.drop_duplicates()
    
    #Detect exons/rows with missing gene symbols or non-numeric values in the percentage30 column
    dataframe['percentage30'] = pd.to_numeric(dataframe['percentage30'], errors='coerce')
    mask_pct30 = dataframe['percentage30'].isna()
    
    if len(dataframe[mask_pct30]) > 0:
        print(f"Warning: {len(dataframe[mask_pct30])} non-numeric values found in 'percentage30' column, removing them from the dataset")
        dataframe = dataframe[~mask_pct30]
    mask_gene = dataframe['GeneSymbol'] == "."

    if len(dataframe[mask_gene]) > 0:
        print(f"Warning: {len(dataframe[mask_gene])} missing gene symbols found in 'GeneSymbol' column, removing them from the dataset")
        dataframe = dataframe[~mask_gene]

    return dataframe


def create_report(input_df: pd.DataFrame, gene_list: str, threshold_value: float) ->  pd.DataFrame:
    """
    Using minimum coverage because it is a more appropriate method to detect short single exon below threshold
    for diagnsotic purposes
    Mean-coverage is good when looking at overall QC metric

    Parameters
    ----------
    input_df : pd.DataFrame
        Pandas dataframe

    gene_list : str
        Name file for the panel_genes.txt

    threshold_value : float
        Threshold value for the coverage30x   

    Returns
    ----------
    pd.DataFrame
        Output report as a pd.DataFrame with columns:
        GeneSymbol  TotalExons  FailingExons  MinCoverage%  WorstExon(chr:start-end)  Status
    """
    print("---Create report---")
    gene_df = pd.read_csv(gene_list, sep=r'\s+', header=None)
    gene_df.columns = ["GeneSymbol"]
    
    results = []
    for gene in gene_df['GeneSymbol']:
    # filter input_df for the current gene
        gene_data = input_df[input_df['GeneSymbol'] == gene]
        if gene_data.empty:
            results.append({
            'GeneSymbol': gene,
            'TotalExons': "NA",
            'FailingExons': "NA",
            'MinCoverage%': "NA",
            'WorstExon(chr:start-end)': "NA",
            'Status': "MISSING"})
        
        else:
            # number of exons per gene
            num_exons = len(gene_data)
    
            # number of failing exons (coverage30x < threshold)
            num_failing = len(gene_data[gene_data['percentage30'] < float(threshold_value)])
    
            # min coverage
            min_coverage = gene_data['percentage30'].min()

            #worst exon:
            min_coverage_idx = gene_data['percentage30'].idxmin()
            worst_exon = gene_data.loc[min_coverage_idx,'FullPosition'].replace('-', ':', 1)

            #status
            status = 'FAIL' if float(min_coverage) < float(threshold_value) else 'PASS'
            results.append({
                'GeneSymbol': gene,
                'TotalExons': num_exons,
                'FailingExons': num_failing,
                'MinCoverage%': min_coverage,
                'WorstExon(chr:start-end)': worst_exon,
                'Status': status})
    
    result_df = pd.DataFrame(results)
    return result_df
    


def save_output(input_df: pd.DataFrame, sambamba_filename:str, threshold_value: float, output_filename):
    """
    Save report file into output.txt and print the run-level summary line
    "Sample NGS148_34: 3/84 genes below 30x threshold"
    Parameters
    ----------
    input_df : pd.DataFrame
        Pandas dataframe

    sambamba_filename: str

    threshold_value: float

    output_filename: str
    """
    print("---Save report---")
    file_name = os.path.basename(sambamba_filename)
    sample_id = '_'.join(file_name.split('_')[0:2])
    mask_failed_gene = (input_df['Status'] == "FAIL")
    num_failed_gene = len(input_df[mask_failed_gene])
    tot_gene = len(input_df)
    print(f"Sample {sample_id}: {num_failed_gene}/{tot_gene} genes below 30x threshold")

    input_df.to_csv(output_filename, index=False, sep="\t")

def exons_failing_gene(sambamba_df, report_df, threshold):
    print("---Write failed_genes_exons.txt with failed exons for failed genes---")
    failed_gene = report_df[(report_df['Status'] == "FAIL")]
    exons = pd.DataFrame()
    for gene in failed_gene['GeneSymbol']:
        gene_exon = sambamba_df[sambamba_df['GeneSymbol'] == gene]
        failed_exons = gene_exon[gene_exon['percentage30'] < threshold]
        exons = pd.concat([exons, failed_exons], ignore_index=True)
    exons.to_csv("failed_genes_exons.txt", index=False, sep="\t")

def main():
    args = parse_args()
    sambamba_df = read_sambamba(args.input)
    output_df = create_report(sambamba_df, args.genelist, args.threshold)
    save_output(output_df, args.input, args.threshold, args.output)
    exons_failing_gene(sambamba_df, output_df, args.threshold)

if __name__ == "__main__":
    main()