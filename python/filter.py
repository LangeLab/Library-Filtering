import os
import re
import time
import argparse
import numpy as np
import pandas as pd

# Example Usages
# --------------
# Minimal: python filter.py library.tsv MsFragger /path/to/main/directory
# Verbose: python filter.py library.tsv MsFragger /path/to/main/directory --verbose
# Output Directory: python filter.py library.tsv MsFragger /path/to/main/directory --output_dir output
# Output File: python filter.py library.tsv MsFragger /path/to/main/directory --output_file filtered_library.tsv
# File Format: python filter.py library.tsv MsFragger /path/to/main/directory --file_format csv
# --------------

# Stores Specific Variables for each library type
library_dict = {
    'Spectronaut': {
        'mod_pattern': '(UniMod:3)',
        'ion_col_name': 'F.FrgIon',
        'charge_col_name': 'FG.Charge',
        'mod_peptide_col_name': 'PEP.GroupingKey',
        'sample_col_name': 'R.FileName',
        'pqColName': 'F.PredictedRelativeIntensity',
        'msColName': 'F.MeasuredRelativeIntensity',
    },
    'MsFragger': {
        'mod_pattern': '(UniMod:3)',
        'ion_col_name': 'Annotation',
        'charge_col_name': 'PrecursorCharge',
        'mod_peptide_col_name': 'ModifiedPeptideSequence',      
        'sample_col_name': None,
        'pqColName': None,
        'msColName': None,
    }
}

# Timer Related Functions
def getTime() -> float:
    """Get the current time for timer."""
    return time.time()

def prettyTimer(seconds: float) -> str:
    """Better way to show elapsed time."""
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%02dh:%02dm:%02ds" % (h, m, s)

def read_spectral_library(
        file_path: str, 
        file_format: str = 'tsv'
    ) -> pd.DataFrame:
    """
    Read a spectral library from a file.

    Args:
        file_path (str): Path to the spectral library file.
        file_format (str): File format ('csv', 'tsv', 'txt', or 'excel').

    Returns:
        pd.DataFrame: A DataFrame containing the loaded spectral library.
    
    Raises:
        ValueError: If an unsupported file format is provided.
    """
    if file_format == 'csv':
        return pd.read_csv(
            file_path, 
            sep=","
        )
    elif file_format == 'tsv' or file_format == 'txt':
        return pd.read_csv(
            file_path, 
            sep="\t"
        )
    elif file_format == 'excel':
        return pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format")

def write_spectral_library(
        library: pd.DataFrame,
        file_path: str,
        file_format: str = 'tsv'
    ) -> None:
    """
    Write the filtered spectral library to a file.
    
    Args:
        library (pd.DataFrame): The filtered spectral library DataFrame.
        file_path (str): Path to the output file.
        file_format (str): File format ('csv', 'tsv', 'txt', or 'excel').

    Raises:
        ValueError: If an unsupported file format is provided.
    """
    if file_format == 'csv':
        library.to_csv(
            file_path,
            sep=",",
            index=False
        )
    elif file_format == 'tsv' or file_format == 'txt':
        library.to_csv(
            file_path,
            sep="\t",
            index=False
        )
    elif file_format == 'excel':
        library.to_excel(
            file_path,
            index=False
        )
    else:
        raise ValueError("Unsupported file format")
        
def process_spectral_library(
        library: pd.DataFrame, 
        library_type: str, 
        verbose: bool = False
    ) -> pd.DataFrame:
    """
    Process the spectral library data.

    Args:
        library (pd.DataFrame): The spectral library DataFrame.
        library_type (str): Type of spectral library ('Spectronaut' or 'MsFragger').
        verbose (bool): Whether to print verbose output (default: False).

    Returns:
        pd.DataFrame: The processed and filtered spectral library DataFrame.
    """

    # Collect the library-specific information from the setup dictionary
    mod_pattern = library_dict[library_type]['mod_pattern']
    ion_col_name = library_dict[library_type]['ion_col_name']
    charge_col_name = library_dict[library_type]['charge_col_name']
    mod_peptide_col_name = library_dict[library_type]['mod_peptide_col_name']
    
    # Following columns are only used for Spectronaut
    if library_type == 'Spectronaut':
        sample_col_name = library_dict[library_type]['sample_col_name']
        pqColName = library_dict[library_type]['pqColName']
        msColName = library_dict[library_type]['msColName']
    

    if verbose:
        print(f"Filtering library from {library_type}...")

    # Initialize the valid_precursors list
    valid_precursors  = []

    # Create a new column to identify precursor
    if library_type == 'Spectronaut':
        library["ID"] = (
            library[mod_peptide_col_name] + "_" + 
            library[sample_col_name] + "_" + 
            library[charge_col_name].astype(str)
        )
        # Specific binary check for Spectronaut
        library["isWithin"] = (
            library[msColName] > (library[pqColName] / 10)
        )
        cols1 = [
            "ID", "isWithin", 
            mod_peptide_col_name, 
            ion_col_name, 
            charge_col_name
        ]
        cols2 = [
            'ID', 'isWithin', 
            'StrippedSequence', 
            'SequenceLength', 
            'ModifiedSequence', 
            'FragmentType', 
            'FragmentNumber'
        ]
    else:
        library["ID"] = (
            library[mod_peptide_col_name] + 
            "_" + 
            library[charge_col_name].astype(str)
        )
        cols1 = [
            "ID", 
            mod_peptide_col_name, 
            ion_col_name, 
            charge_col_name
        ]
        cols2 = [
            'ID', 
            'StrippedSequence', 
            'SequenceLength', 
            'ModifiedSequence', 
            'FragmentType', 
            'FragmentNumber'
        ]

    if verbose:
        # Print the number of unique precursors in the initial data
        print(f"Initial Data - Number of Unique Precursors: {len(set(library['ID']))}")

    # Select a subset of the data
    data = library[cols1].copy()

    # N-term Modification pattern 
    #   Matches the beginning of the string and 
    #   Matches the pattern of the form `(UniMod:[0-9]+)`
    # TODO: In the future different labeling styles need to be considered
    nterm_mod_pattern = re.compile(r'^\([^)]*\)') 
    # Remove all rows that do not have a modification at the N-terminus
    data = data[
        data[mod_peptide_col_name].astype(str).apply(
            lambda x: bool(nterm_mod_pattern.search(x))
        )
    ]

    if verbose:
        # Print the number of unique precursors after selecting Nterm label only precusors
        print(f"Labeled at Nterm Only - Number of Unique Precursors: {len(set(data['ID']))}")

    # Get Fragment Type and Number
    data[[
        'FragmentType',
        'FragmentNumber'
    ]] = data[ion_col_name].str.split(
        "^" # Removes the trailing ^ and the number after it
    ).str[0].str.extract(
        r'([a-zA-Z]+)(\d+)' # Extracts fragment and number
    )
    # Convert FragmentNumber to int
    data["FragmentNumber"] = data["FragmentNumber"].astype(int)

    # Create new columns to store modified (without N-term mod) and stripped sequences
    data["ModifiedSequence"] = data[mod_peptide_col_name].str.replace(
        # Only remove the N-term mod
        nterm_mod_pattern, 
        ''
    )
    data["StrippedSequence"] = data["ModifiedSequence"].str.replace(
        # Removes all the mod labels with the form `(UniMod:[0-9]+)`
        r'\(UniMod:\d+\)', 
        '', 
        regex=True
    )
    data["SequenceLength"] = data["StrippedSequence"].str.len()

    # Only keep the columns we need
    data = data[cols2]

    # Find the position of all Ks on the stripped sequence
    data["KPositions"] = data.apply(
        lambda row: [
            m.start() for m in re.finditer(
                r'K', row["StrippedSequence"]
            )
        ], 
        axis=1
    )

    # Subset for all rows that do not have any K
    subset = data[
        data["KPositions"].str.len() == 0
    ].reset_index(drop=True).copy()

    # Place the valid precursors without K in the valid_precursors list
    if library_type == 'Spectronaut':
        valid_precursors.extend(
            subset[
                subset["isWithin"]
            ]["ID"].unique()
        )
    else:
        valid_precursors.extend(
            subset["ID"].unique()
        )
    if verbose:
        print(f"Number of Unique Valid Precursors (without any K): {len(valid_precursors)}")

    # Remove all rows without any K since we already checked them
    data = data[
        data["KPositions"].str.len() > 0
    ].reset_index(drop=True).copy() 

    # Find all the modification on the modified sequence
    data["Modifications"] = data["ModifiedSequence"].str.findall(
        r'\(UniMod:[0-9]+\)' # TODO: More robust to multiple mod labeling types
    )
    # Find the positions of the modifications
    data["Positions"] = data.apply(
        lambda row: [
            m.start() - 1 for m in re.finditer(
                r'\(UniMod:[0-9]+\)',  # TODO: More robust to multiple mod labeling types
                row["ModifiedSequence"]
            )
        ], 
        axis=1
    )
    # Find the preceding amino acids per modification
    data["PrecedingAminoAcids"] = data.apply(
        lambda row: [
            row["ModifiedSequence"][pos] if pos >= 0 else None for pos in row["Positions"]
        ], 
        axis=1
    )
    # Find the absolute (without mod labels) positions
    data["AbsolutePositions"] = data.apply(
        lambda row: [
            pos - cumulative_length for pos, cumulative_length in zip(
                row["Positions"], 
                (
                    np.cumsum([len(mod) for mod in row["Modifications"]]) - 
                    [len(mod) for mod in row["Modifications"]]
                )
            )
        ],
        axis=1
    )

    # Find the K positions after the fragment number for y-ion check
    data["AfterFragmentKPos"] = data.apply(
        lambda row: [
            p for p in row["KPositions"] if p >= (row["SequenceLength"] - int(row["FragmentNumber"]))
        ], 
        axis=1
    )
    # Find the K positions before the fragment number for b-ion check
    data["BeforeFragmentKPos"] = data.apply(
        lambda row: [
            p for p in row["KPositions"] if p < int(row["FragmentNumber"]) 
        ], 
        axis=1
    )
    # Find the modified K positions using the absolute positions
    data["ModifiedKPositions"] = data.apply(
        lambda row: [
            p for p in row["KPositions"] if p in row["AbsolutePositions"]
        ], 
        axis=1
    )
    # Find the unmodified K positions using the absolute positions
    data["UnmodifiedKPositions"] = data.apply(
        lambda row: [
            p for p in row["KPositions"] if p not in row["AbsolutePositions"]
        ], 
        axis=1
    )
    # Check if the fragment ion and the number are correct for K and Biotin Labels
    data["correct_fragment_ion"] = data.apply(
        lambda row: (
            (
                row["FragmentType"] == 'b' and 
                set(row["BeforeFragmentKPos"]).issubset(set(row["ModifiedKPositions"]))
            ) or
            (
                row["FragmentType"] == 'y' and 
                set(row["UnmodifiedKPositions"]).issubset(set(row["AfterFragmentKPos"]))
            )
        ), 
        axis=1
    )

    # Extend the valid_precursors list with the valid precursors corresponding to correct_fragment_ion
    if library_type == 'Spectronaut':
        valid_precursors.extend(
            data[
                data["correct_fragment_ion"] & 
                data["isWithin"]
            ]["ID"].unique()
        )
    else:
        valid_precursors.extend(
            data[
                data["correct_fragment_ion"]
            ]["ID"].unique()
        )
    if verbose:
        print(f"Number of Unique Valid Precursors (all): {len(valid_precursors)}")

    # Filter the library using the valid_precursors list
    filtered_library = library[
        library["ID"].isin(valid_precursors)
    ].drop(columns=["ID"])

    # Return the filtered library
    return filtered_library

def main():
    """
    Main function for processing spectral library data.
    
    Parameters:
        None
        
    Returns:
        None
    """
    # Start timer
    start_time = getTime()
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process spectral library data")
    parser.add_argument("file_name", help="Name of the spectral library file")
    parser.add_argument("library_type", choices=["Spectronaut", "MsFragger"], help="Type of spectral library")
    parser.add_argument("main_path", help="Main directory where the input file is located")
    parser.add_argument("--output_dir", help="Output directory for the filtered library")
    parser.add_argument("--output_file", help="Specify the output file name for the filtered library")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    args = parser.parse_args()

    # Get the variables from the arguments
    file_name = args.file_name
    library_type = args.library_type
    main_path = args.main_path
    output_dir = args.output_dir
    output_file = args.output_file
    verbose = args.verbose

    # Find the file format from the file name
    file_format = file_name.split(".")[-1]
    file_root_name = ".".join(file_name.split(".")[:-1])
    
    # Create check for library type
    if library_type not in ['Spectronaut', 'MsFragger']:
        raise ValueError("Unsupported library type, please choose from Spectronaut or MsFragger")
    # Create check for file path
    if not os.path.exists(os.path.join(main_path, file_name)):
        raise ValueError("File is not found in the specified path with filename")
    if verbose:
        print(f"Reading {file_name} from {main_path}...")

    # Construct the full file path using main_path
    full_file_path = os.path.join(main_path, file_name)

    # Read the spectral library
    library = read_spectral_library(full_file_path, file_format)
    # Apply the filter to the library
    filtered_library = process_spectral_library(
        library,
        library_type,
        verbose
    )

    # If output_file is not specified, use the input file name
    if not output_file:
        # update the file_name with _filtered (more robust way to do this?)
        output_file = file_root_name + "_filtered." + file_format

     # If output_dir is specified, use it as the base directory
    if output_dir:
        output_dir = os.path.join(main_path, output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_file_path = os.path.join(output_dir, output_file)
    else:
        # Use main_path for output
        output_file_path = os.path.join(main_path, output_file)

    # Save the filtered library
    write_spectral_library(filtered_library, output_file_path, file_format)

    elapsed_time = getTime() - start_time
    print(f"Filtering complete! Elapsed Time: {prettyTimer(elapsed_time)}")

if __name__ == "__main__":
    main()