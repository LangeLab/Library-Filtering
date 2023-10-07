suppressMessages(library(dplyr))  # for data manipulation
suppressMessages(library(readr))  # for reading csv files
suppressMessages(library(readxl))  # for reading excel files
suppressMessages(library(writexl))  # for writing excel files
suppressMessages(library(stringr))  # for string manipulation
suppressMessages(library(optparse))  # for command line arguments

# Example Usage:
# Rscript filter.R -f MSFragger_spectral_library.tsv -s MsFragger -m ../Data/test_data/
# Rscript filter.R -f Spectronaut_report.tsv -s Spectronaut -m ../Data/test_data/ 

debug_mode <- FALSE

# Stores Specific Variables for each library type
library_dict <- list(
  Spectronaut = list(
    mod_pattern = "(UniMod:3)",
    ion_col_name = "F.FrgIon",
    charge_col_name = "FG.Charge",
    mod_peptide_col_name = "PEP.GroupingKey",
    sample_col_name = "R.FileName",
    pqColName = "F.PredictedRelativeIntensity",
    msColName = "F.MeasuredRelativeIntensity"
  ),
  MsFragger = list(
    mod_pattern = "(UniMod:3)",
    ion_col_name = "Annotation",
    charge_col_name = "PrecursorCharge",
    mod_peptide_col_name = "ModifiedPeptideSequence",
    sample_col_name = NULL,
    pqColName = NULL,
    msColName = NULL
  )
)

getTime <- function() {
  return(Sys.time())
}

prettyTimer <- function(seconds) {
  seconds <- as.numeric(seconds)
  m <- floor(seconds / 60)
  s <- seconds %% 60
  h <- floor(m / 60)
  m <- m %% 60
  return(sprintf("%02.0fh:%02.0fm:%02.0fs", h, m, s))
}

read_spectral_library <- function(file_path, file_format = "tsv") {
  if (file_format == "csv") {
    return(readr::read_csv(file_path, show_col_types = FALSE))
  } else if (file_format %in% c("tsv", "txt")) {
    return(readr::read_tsv(file_path, show_col_types = FALSE))
  } else if (file_format == "excel") {
    return(readxl::read_excel(file_path))
  } else {
    stop("Unsupported file format")
  }
}

write_spectral_library <- function(library, file_path, file_format = "tsv") {
  if (file_format == "csv") {
    readr::write_csv(library, file_path)
  } else if (file_format %in% c("tsv", "txt")) {
    readr::write_tsv(library, file_path)
  } else if (file_format == "excel") {
    writexl::write_xlsx(library, file_path)
  } else {
    stop("Unsupported file format")
  }
}

process_spectral_library <- function(library, library_type, verbose = FALSE) {
  
  # Collect the library-specific information from the setup dictionary
  mod_pattern <- library_dict[[library_type]]$mod_pattern
  ion_col_name <- library_dict[[library_type]]$ion_col_name
  charge_col_name <- library_dict[[library_type]]$charge_col_name
  mod_peptide_col_name <- library_dict[[library_type]]$mod_peptide_col_name
  
  # Following columns are only used for Spectronaut
  if (library_type == 'Spectronaut') {
    sample_col_name <- library_dict[[library_type]]$sample_col_name
    pqColName <- library_dict[[library_type]]$pqColName
    msColName <- library_dict[[library_type]]$msColName
  }
  
  if (verbose) {
    message(paste("Filtering the file (source:", library_type, ")"))
  }

  # Initialize the valid_precursors list
  valid_precursors <- c()
  
  # Create a new column to identify precursor
  if (library_type == 'Spectronaut') {
    # Create identifier
    library$ID <- paste0(
      library[[mod_peptide_col_name]], "_",
      library[[sample_col_name]], "_",
      as.character(library[[charge_col_name]])
    )
    # Create a binary checker
    library$isWithin <- library[[msColName]] > (library[[pqColName]] / 10)
    # Specify subsetter vectors
    cols1 = c('ID', 'isWithin', mod_peptide_col_name, ion_col_name, charge_col_name)
    cols2 = c('ID', 'isWithin', 'StrippedSequence', 'SequenceLength', 'ModifiedSequence', 'FragmentType', 'FragmentNumber')
  } else {
    # Create identifier
    library$ID <- paste0(
      library[[mod_peptide_col_name]], "_",
      as.character(library[[charge_col_name]])
    )
    cols1 = c("ID", mod_peptide_col_name, ion_col_name, charge_col_name)
    cols2 = c('ID', 'StrippedSequence', 'SequenceLength', 'ModifiedSequence', 'FragmentType', 'FragmentNumber')
  }

  if (debug_mode){
    message(paste("Full library has", nrow(library), "rows"))
  }

  if (verbose) {
    # Print the number of unique precursors in the initial data
    message(paste(" -> Initial Data has", length(unique(library$ID)), "unique ids to check"))
  }
  
  # N-term Modification pattern
  # Matches the beginning of the string and
  # Matches the pattern of the form `(UniMod:[0-9]+)`
  # TODO: In the future different labeling styles need to be considered
  nterm_mod_pattern <- "^\\([^)]*\\)"

  # Select a subset of the data
  data <- library %>% 
    # Select only columns to be used
    select(all_of(cols1)) %>%
    # Filter entries that start with nterm_mod_pattern
    filter(str_detect(.data[[mod_peptide_col_name]], nterm_mod_pattern))
  
  if (verbose) {
    # Print the number of unique precursors after selecting Nterm label only precursors
    message(paste(" -> Labeled at Nterm Only has", length(unique(data$ID)), "unique ids to check"))
  }
  
  # Get Fragment Type and Number
  data <- data %>%
    mutate(
      FragmentType = substr(.data[[ion_col_name]], 1, 1),
      FragmentNumber = as.integer(str_extract(.data[[ion_col_name]], "\\d+"))
    )  

  # Create new columns to store modified (without N-term mod) and stripped sequences
  data <- data %>%
    mutate(
      # Remove the first mod_pattern (which all rows start with)
      ModifiedSequence = str_remove(.data[[mod_peptide_col_name]], nterm_mod_pattern),
      StrippedSequence = str_replace_all(
        ModifiedSequence,
        "\\(UniMod:\\d+\\)",
        ''
      ),
      SequenceLength = str_length(StrippedSequence)
    ) %>% 
    select(all_of(cols2))

  # Find the position of all Ks on the stripped sequence
  data$KPositions <- sapply(
    data$StrippedSequence, 
    function(x) which(strsplit(x, "")[[1]] == "K")
  )

  # If KPositions is empty, replace it with NA
  data$KPositions[data$KPositions == "integer(0)"] <- NA
    
  if (debug_mode) {
    message(paste("N-term Only data has", nrow(data), "rows"))
  }
  # Subset the data if there are no Ks
  subset <- data[is.na(data$KPositions), ]
    
  if (debug_mode) {
    message(paste("N-term Only & No-K Sequence subset has", nrow(subset), "rows"))
  }
  

# Subset for no Ks only
subset <- data[is.na(data$KPositions), ]

# If only n-term labeled without any Ks in the sequence all are valid
if (library_type == 'Spectronaut') {
  valid_precursors <- unique(
    c(
      valid_precursors, 
      subset$ID[subset$isWithin]
    )
  )
} else {
  valid_precursors <- unique(
    c(
      valid_precursors, 
      subset$ID
    )
  )
}
  if (verbose) {
    # Print the number of unique precursors after No-K filtering
    message(paste(" -> After checking No-K sequences", length(valid_precursors), "valid unique ids found"))
  }

  # Select a subset of the data that hase at least one K residue
  data <- data[!is.na(data$KPositions), ]

  if (debug_mode) {
    message(paste("K-Only Sequence subset has", nrow(data), "rows"))
  }

  # Find all the modification on the modified sequence
  data$Modifications <- regmatches(
    data$ModifiedSequence, 
    gregexpr(
      "\\(UniMod:[0-9]+\\)", 
      data$ModifiedSequence
    )
  )
  # If Modifications is empty, replace it with NA
  data$Modifications[data$Modifications == "character(0)"] <- NA

  # Find the positions of the modifications
  data$Positions <- mapply(
    function(seq) {
      match <- regmatches(seq, gregexpr("\\(UniMod:[0-9]+\\)", seq))
      if (length(match[[1]]) == 0) {
        NA
      } else {
        as.integer(gregexpr("\\(UniMod:[0-9]+\\)", seq)[[1]] - 1)
      }
    },
    data$ModifiedSequence
  )

  # Find the preceding amino acids per modification
  data$PrecedingAminoAcids <- mapply(
    function(seq, positions) {
      sapply(
        positions, 
        function(pos) ifelse(pos >= 0, substr(seq, pos, pos), NA)
      )
    }, 
    data$ModifiedSequence, 
    data$Positions
  )

  # Find the absolute (without mod labels) positions
  data$AbsolutePositions <- mapply(
    function(positions, mods) {
      positions - cumsum(sapply(mods, nchar)) + sapply(mods, nchar)
    }, 
    data$Positions, 
    data$Modifications
  )

  # Find the K positions after the fragment number for y-ion check
  data$AfterFragmentKPos <- mapply(
    function(k_positions, fragment_number, seq_length) {
      k_positions[k_positions > (seq_length - fragment_number)]
    }, 
    data$KPositions, 
    data$FragmentNumber, 
    data$SequenceLength
  )
  # If AfterFragmentKPos is empty, replace it with NA
  data$AfterFragmentKPos[data$AfterFragmentKPos == "integer(0)"] <- NA

  # Find the K positions before the fragment number for b-ion check
  data$BeforeFragmentKPos <- mapply(
    function(k_positions, fragment_number) {
      k_positions[k_positions <= (fragment_number)]
    }, 
    data$KPositions, 
    data$FragmentNumber
  )
  # If BeforeFragmentKPos is empty, replace it with NA
  data$BeforeFragmentKPos[data$BeforeFragmentKPos == "integer(0)"] <- NA

  # Find the modified K positions using the absolute positions
  data$ModifiedKPositions <- mapply(
    function(k_positions, absolute_positions) {
      k_positions[k_positions %in% absolute_positions]
    }, 
    data$KPositions, 
    data$AbsolutePositions
  )
  # If ModifiedKPositions is empty, replace it with NA
  data$ModifiedKPositions[data$ModifiedKPositions == "integer(0)"] <- NA

  # Find the unmodified K positions using the absolute positions
  data$UnmodifiedKPositions <- mapply(
    function(k_positions, absolute_positions) {
      k_positions[!(k_positions %in% absolute_positions)]
    }, 
    data$KPositions, 
    data$AbsolutePositions
  )
  # If UnmodifiedKPositions is empty, replace it with NA
  data$UnmodifiedKPositions[data$UnmodifiedKPositions == "integer(0)"] <- NA

  # Check if the fragment ion and the number are correct for K and Biotin Labels
  data$correct_fragment_ion <- mapply(
    function(fragment_type, before_k, after_k, unmodified_k, modified_k) {
      (fragment_type == 'b' && (all(is.na(before_k)) || all(before_k %in% modified_k))) ||
      (fragment_type == 'y' && all(unmodified_k %in% after_k))
    }, 
    data$FragmentType, 
    data$BeforeFragmentKPos, 
    data$AfterFragmentKPos, 
    data$UnmodifiedKPositions,
    data$ModifiedKPositions
  ) 

  # Extend the valid_precursors list with the valid precursors corresponding to correct_fragment_ion
  if (library_type == 'Spectronaut') {
    valid_precursors <- unique(
      c(
        valid_precursors, 
        data$ID[data$correct_fragment_ion & data$isWithin]
      )
    )
  } else {
    valid_precursors <- unique(
      c(
        valid_precursors, 
        data$ID[data$correct_fragment_ion]
      )
    ) 
  }

  if (verbose) {
    # Print the number of unique precursors after processing
    message(paste(" -> After filtering there are", length(unique(data$ID)), "unique ids are valid"))
  }

  # Filter the library using the valid_precursors list
  filtered_library <- library[library$ID %in% valid_precursors, , drop = FALSE]
  filtered_library <- filtered_library[, !names(filtered_library) %in% "ID", drop = FALSE]
  
  return(filtered_library)
}
  
main <- function() {
  # Start timer
  start_time <- getTime()

  # Define the options
  option_list <- list(
    make_option(
      c("-f", "--file"), 
      type = "character", 
      help = "Name of the spectral library file",
    ),
    make_option(
      c("-s", "--lib_source"), 
      type = "character", 
      help = "Type of spectral library",
    ),
    make_option(
      c("-m", "--main_path"), 
      type = "character", 
      help = "Main directory where the input file is located",
    ),
    make_option(
      c("-o", "--output_dir"), 
      type = "character", 
      help = "Output directory for the filtered library",
      default = NULL
    ),
    make_option(
      c("-r", "--output_file"), 
      type = "character", 
      help = "Specify the output file name for the filtered library",
      default = NULL
    ),
    make_option(
      c("-v", "--verbose"), 
      action = "store_true", 
      help = "Print verbose output", 
      default = FALSE
    )
  )

  # Parse the command-line arguments
  opt_parser <- OptionParser(usage = "Rscript filter.R [options]", option_list = option_list)
  opt <- parse_args(opt_parser)

  # Get the variables from the options
  file_name <- opt$file
  library_type <- opt$lib_source
  main_path <- opt$main_path
  output_dir <- opt$output_dir
  output_file <- opt$output_file
  verbose <- opt$verbose

  # Check if 'type' is valid
  if (!(library_type %in% c("Spectronaut", "MsFragger"))) {
      stop("Error: The --type option must be either Spectronaut or MsFragger.")
  }

  # Find the file format from the file name
  file_format <- strsplit(file_name, "\\.")[[1]][length(strsplit(file_name, "\\.")[[1]])]
  file_root_name <- paste(strsplit(file_name, "\\.")[[1]][-length(strsplit(file_name, "\\.")[[1]])], collapse=".")

  # Create check for file path
  if (!file.exists(file.path(main_path, file_name))) {
      stop("File is not found in the specified path with filename")
  }

  if (verbose) {
      message(paste0("Reading ", file_name, " from ", main_path))
  }

  # Construct the full file path using main_path
  full_file_path <- file.path(main_path, file_name)

  # Read the spectral library
  library <- read_spectral_library(full_file_path, file_format)
  
  # Apply the filter to the library
  filtered_library <- process_spectral_library(
      library,
      library_type,
      verbose
  )

  # If output_file is not specified, use the input file name
  if (is.null(output_file)) {
      # update the file_name with _filtered (more robust way to do this?)
      output_file <- paste(file_root_name, "_filtered.", file_format, sep="")
  }

  # If output_dir is specified, use it as the base directory
  if (!is.null(output_dir)) {
      output_dir <- file.path(main_path, output_dir)
      if (!dir.exists(output_dir)) {
          dir.create(output_dir)
      }
      output_file_path <- file.path(output_dir, output_file)
  } else {
      # Use main_path for output
      output_file_path <- file.path(main_path, output_file)
  }

  # Save the filtered library
  write_spectral_library(filtered_library, output_file_path, file_format)

  elapsed_time <- getTime() - start_time
  message(paste0("Filtering complete! Elapsed Time: ", prettyTimer(elapsed_time)))
}

# Run the main function
main()
