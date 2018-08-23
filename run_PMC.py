from peptide_mass_calculator import *


# Tries to open seuqnece file and build the required objects for script to round
# Creates a sequence file otherwise
while True:
    try:
        aaIonMasses = getAminoAcidMasses()
        print('Loading data from json...')
        sequence = getSequenceAAs()
        print('Sequence file found, reading now...')
        uniqueAA = getUniqueAA(sequence)
        print('Getting unique amino acids...')
        mass = getExactMass(sequence, aaIonMasses)
        print('Calculating mass...')
        possibleMasses = populateMassDictionary(mass, sequence, uniqueAA, aaIonMasses)
        print('Calculating possible ions...')
        break
    except IOError:
        print("No sequence file found, creating now...")
        createSequenceFile()


# User searches for a match, ensures that the user enters a float
while True:
    try:
        matches = getMassMatch(float(input("Enter your observed mass: ")), possibleMasses)
        for match in matches:
            print(match)
        break
    except ValueError:
        print("please enter your mass in the form ####.#### or ####")

# Allows user to enter another match, or enter
while True:
    try:
        matches = getMassMatch(float(input("Enter another observed mass, or hit enter key to exit: ")), possibleMasses)
        for match in matches:
            print(match)
    except:
        writeMassDictionary(possibleMasses)
        quit()
