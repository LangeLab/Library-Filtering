import argparse
import pandas as pd
import re

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
    modPeptideColName = 'ModifiedPeptide'
    stripPeptideColName = 'StrippedPeptide'
    ionNumColName = 'FragmentNumber'
elif libraryType == 'DiaNN':
    # filter out rows w/ no (Unimod:3) at the start (N-terminus) of ModifiedPeptideSequence
    library = library[library.apply(lambda row: row['ModifiedPeptideSequence'][:10] == '(UniMod:3)', axis=1)]
    modPeptideColName = 'ModifiedPeptideSequence'
    stripPeptideColName = 'PeptideSequence'
    ionNumColName = 'FragmentSeriesNumber'
else:
    print('Library type must be "Spectronaut" or "DiaNN".')
    exit()

# nested dictionary that will store for each peptide:
# - list of the 0-based positions of the modified (biotinylated) K residues
# - list of the 0-based positions of the unmodified (not biotinylated) K residues
# - first confirming ion detected for a peptide of distinct modified sequence and charge
modifiedPeptides = {}

def getModKPositions(row):
    if libraryType == 'DiaNN':
        # get a version of the ModifiedPeptideSequence that has K(UniMod:3) as the only modification
        string = re.sub('[^K]\(UniMod:[0-9]+\)', 'X', row['ModifiedPeptideSequence'][10:])
        string = re.sub('K\(UniMod:([^3]|[0-9]{2,})\)', 'K', string)
        # find all 0-based positions of K(UniMod:3)
        indices = [index for index in range(len(string)) if string.startswith('K(UniMod:3)', index)]
        # take into account that (UniMod:3) substrings are shifting index values
        indices = [value - (10 * index) for index, value in enumerate(indices)]
        return indices
    else:
        # library type is Spectronaut
        # get a version of ModifiedPeptide that has [Biotin-NHS-K] as the only modification
        string = re.sub('[^K]\[[^\]]*\]', 'X', row['ModifiedPeptide'].strip('_')[20: ])
        # assuming ? symbol is not originally used anywhere in the ModifiedPeptide notation
        string = string.replace('Biotin-NHS-K', '?')
        string = re.sub('\[[^\?]*\]', '', string)
        string = string.replace('?', 'Biotin-NHS-K')
        indices = [index for index in range(len(string)) if string.startswith('K[Biotin-NHS-K]', index)]
        # take into account that [Biotin-NHS-K] substrings are shifting index values
        indices = [value - (14 * index) for index, value in enumerate(indices)]
        return indices

def check(row):
    keyName = row[modPeptideColName] + '|' + str(row['PrecursorCharge'])
    if keyName not in list(modifiedPeptides.keys()):
        # initialize entry into dict
        modifiedPeptides[keyName] = {'modified_K_positions': getModKPositions(row)}
        modifiedPeptides[keyName]['unmodified_K_positions'] = list(set([i for i in range(len(row[stripPeptideColName])) if row[stripPeptideColName].startswith('K', i)]) - set(modifiedPeptides[keyName]['modified_K_positions'])) 
        modifiedPeptides[keyName]['confirming_fragment_ion'] = ''
    if modifiedPeptides[keyName]['confirming_fragment_ion'] == '':
        if row['FragmentType'] == 'b':
            b_ion_frag = row[stripPeptideColName][: row[ionNumColName]]
            indices = [i for i in range(len(b_ion_frag)) if b_ion_frag.startswith('K', i)]
            # b-ion fragment can have K residues as long as they're modified
            if set(indices).issubset(set(modifiedPeptides[keyName]['modified_K_positions'])): 
                modifiedPeptides[keyName]['confirming_fragment_ion'] = row['FragmentType'] + str(row[ionNumColName])          
        elif row['FragmentType'] == 'y':
            y_ion_frag = row[stripPeptideColName][-row[ionNumColName]: ]
            indices = [len(row[stripPeptideColName]) - len(y_ion_frag) + i for i in range(len(y_ion_frag)) if y_ion_frag.startswith('K', i)]
            # y-ion fragment needs to have at least all of the unmodified K residues           
            if set(modifiedPeptides[keyName]['unmodified_K_positions']).issubset(set(indices)):
                modifiedPeptides[keyName]['confirming_fragment_ion'] = row['FragmentType'] + str(row[ionNumColName])
    
library.apply(lambda row: check(row), axis=1)

# list of peptides to keep
filteredPeptides = []
for peptide in modifiedPeptides:
    if modifiedPeptides[peptide]['confirming_fragment_ion'] != '':
        filteredPeptides.append(peptide)

# filtered library below
library = library[library.apply(lambda row: row[modPeptideColName] + '|' + str(row['PrecursorCharge']) in filteredPeptides, axis=1)]

# output filtered library in the same directory as script
library.to_csv('filtered_' + fileName, sep='\t', index=False)
print('Filtering complete.')
