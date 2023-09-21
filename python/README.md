# Python Implementation of Filtering Algorithm

## Usage

User must specify the input file name, source (determines the quant or qual), input path. The output directory and file name are optional. If not specified, the output directory will be the same as the input path and the output file name will be the same as the input file name with the suffix "_filtered".

```bash
usage: filter.py [-h] [--output_dir OUTPUT_DIR] [--output_file OUTPUT_FILE] [--verbose] file_name {Spectronaut,MsFragger} main_path

Process spectral library data

positional arguments:
  file_name             Name of the spectral library file
  {Spectronaut,MsFragger}
                        Type of spectral library
  main_path             Main directory where the input file is located

optional arguments:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Output directory for the filtered library
  --output_file OUTPUT_FILE
                        Specify the output file name for the filtered library
  --verbose             Print verbose output
```

## Example

```bash
python filter.py data.tsv Spectronaut /path/to/data --output_dir /path/to/output --output_file filtered_data.tsv --verbose
```
