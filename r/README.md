# R Implementation of Filtering Algorithm

## Usage

User must specify the input file name, source (determines the quant or qual), input path. The output directory and file name are optional. If not specified, the output directory will be the same as the input path and the output file name will be the same as the input file name with the suffix "_filtered".

```bash
Usage: Rscript filter.R [options]

Options:
        -f FILE, --file=FILE
                Name of the spectral library file

        -s LIB_SOURCE, --lib_source=LIB_SOURCE
                Type of spectral library

        -m MAIN_PATH, --main_path=MAIN_PATH
                Main directory where the input file is located

        -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                Output directory for the filtered library

        -r OUTPUT_FILE, --output_file=OUTPUT_FILE
                Specify the output file name for the filtered library

        -v, --verbose
                Print verbose output

        -h, --help
                Show this help message and exit
```

## Example

```bash
Rscript filter.R -f data.tsv -s Spectronaut -m /path/to/data -o /path/to/output -r filtered_data.tsv -v
```
