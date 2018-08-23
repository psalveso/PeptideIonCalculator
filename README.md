# PeptideIonCalculator

This python script lets users calculate peptide ions (m/z) that will be observed in ESI-MS. It allows users to enter an observed ion, and tells the user what that ion is (e.g., [M+2H]2+)

## Dependencies

- Python3

## Running the script

After downloading the repository and ensuring the required Dependencies are installed, navigate to the folder containing the source files using the terminal.

To run the script type: `python3 run_PMC.py`

### Building you Sequence
Once the script is running, it will first look for a `sequence.txt` file in the current directory. If it is not there, it will prompt you to make one. Type you amino acid sequence, using the corresponding codes, which are viewable in `aaIonMasses.json`.

Add the following residues to your sequence if your peptide is:

- `HO` for linear peptides
- `Cyclic` for cyclic peptides (lactams or lactones)
- `NH2` for linear peptides with *C*-terminal amides
- `NHCH3` for linear peptides with *C*-terminal methylamides

The script will create the `sequence.txt` file, and can be rerun so that you do not have to constantly re-input the sequence.

### Searching for matches
Once your sequence is entered, the script will prompt you for an observed mass. Enter your mass, and the script will return any matching ion that is within 0.3 amu of your queried mass. For example, if you search for 736, the script will return all ions whose m/z is between 735.7 and 736.3. If you do not get a match, try adjusting your queried mass by 0.5 amu, to accommodate fo a possible miss-calibrated instrument.

#### Ions it will search against
Currently, the script calculates ions of the form: `[l M + i H + j Na + k K]z+`, where `l`, `i`, `j`, `k`, and `z` are integers. For example: `[1M + 2H + K]3+`. It will calculate ions with z no larger than 3 (i.e., triply charged ions) and for aggregates no larger than 2.

In the future, I will add the ability for users to specify their own limits to `l`, `i`, `j`, `k`, and `z`

The script will also calculate a number of side reactions that can occur, based on the sequence of the peptide. For example, if your sequence contains lysine residues, the script will also calculate ions that contain the mass for an incomplete deprotection of the Boc protecting group. It will do this for tBu protecting groups, Acm protecting groups, Trityl protecting groups, as well as possible TFA esters, assuming those modifications or possible given the amino acids that comprise your peptide.

When the script is complete, it will create a `possible matches.txt` file, which contains all of the ions that it has calculated for you sequence.

In the future, I will add single and double deletions, as well as single and double additions to the possible ions you might see in the mass spectrum. These side reactions can help diagnose failed syntheses.

## Adding your own amino acid
Amino acid codes, and their corresponding masses, are stored in `aaIonMasses.json`. Each entry in this `.json` is formatted as a sting:float, where the string is the amino acid (I've used short names) and the float is the mass of the amino acid. The masses in this `.json` are not the molecular weight of the amino acids, but instead are the mass that each amino acid adds to the peptide. In other words, these masses are calculated as: **exact mass of the amino acid - H2O**.

The amino acid codes can be anything, just ensure that no code is repeated for any amino acid. In the future I will add a tool to ensure amino acids are added to  the `.json` properly.
