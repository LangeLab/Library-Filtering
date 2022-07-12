# Library-Filtering
 Script to filter spectral libraries from Spectronaut or DiaNN for peptides with fragments that support biotin modification at the peptide N terminus.

## Usage:
A user must provide the library file name/path and indicate whether the file is from Spectronaut or DiaNN.

```
usage: library_filter.py [-h] [--inputFile INPUTFILE] [--libraryType LIBRARYTYPE]

options:
  -h, --help            show this help message and exit
  --inputFile INPUTFILE
                        name/path of tab-separated file
  --libraryType LIBRARYTYPE
                        one of the following: Spectronaut, DiaNN
```

## Example:
```
py library_filter.py --inputFile library.tsv --libraryType DiaNN
```