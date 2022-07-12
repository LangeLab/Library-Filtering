import argparse
import pandas as pd

# get file name from command line using the --inputFile flag
# get library type from command line using the --libraryType flag
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--inputFile", type=str, help="name/path of tab-separated file")
parser.add_argument("--libraryType", type=str, help="one of the following: Spectronaut, DiaNN")
args = vars(parser.parse_args())

fileName = str(args['inputFile'])
libraryType = str(args['libraryType'])

library = pd.read_csv(fileName, sep='\t')

if libraryType == 'Spectronaut':
    # filter out rows w/ no [Biotin -NHS-N term] in ModifiedPeptide
    library = library[library.apply(lambda row: '[Biotin -NHS-N term]' in row['ModifiedPeptide'], axis=1)]
    peptideColName = 'StrippedPeptide'
    ionNumColName = 'FragmentNumber'
elif libraryType == 'DiaNN':
    # filter out rows w/ no (Unimod:3) at the start (N-terminus) of ModifiedPeptideSequence
    library = library[library.apply(lambda row: row['ModifiedPeptideSequence'][:10] == '(UniMod:3)', axis=1)]
    peptideColName = 'PeptideSequence'
    ionNumColName = 'FragmentSeriesNumber'
else:
    print('Library type must be "Spectronaut" or "DiaNN".')
    exit()

# nested dictionary that will store for each peptide:
# - the 1-based position of the first K (or 0 if no lysines present)
# - total number of lysines (K) in peptide
# - whether it meets the criteria to keep
modifiedPeptides = {}

def check(row):
    if row[peptideColName] not in list(modifiedPeptides.keys()):
        # initialize entry into dict
        modifiedPeptides[row[peptideColName]] = {'first_K_pos': row[peptideColName].find('K') + 1}
        modifiedPeptides[row[peptideColName]]['total_num_K'] = row[peptideColName].count('K')
        modifiedPeptides[row[peptideColName]]['meetsCriteria'] = False
    if not modifiedPeptides[row[peptideColName]]['meetsCriteria']:
        if row['FragmentType'] == 'b':
            b_ion_frag = row[peptideColName][: row[ionNumColName]]
            if 'K' not in b_ion_frag:
                modifiedPeptides[row[peptideColName]]['meetsCriteria'] = True
        elif row['FragmentType'] == 'y':
            y_ion_frag = row[peptideColName][-row[ionNumColName]: ]
            if y_ion_frag.count('K') == modifiedPeptides[row[peptideColName]]['total_num_K']:
                modifiedPeptides[row[peptideColName]]['meetsCriteria'] = True
    
library.apply(lambda row: check(row), axis=1)

# list of peptides to keep
filteredPeptides = []
for peptide in modifiedPeptides:
    if modifiedPeptides[peptide]['meetsCriteria']:
        filteredPeptides.append(peptide)

# filtered library below
library = library[library.apply(lambda row: row[peptideColName] in filteredPeptides, axis=1)]

# output filtered library in the same directory as script
library.to_csv('filtered_' + fileName, sep='\t', index=False)
print('Filtering complete.')
