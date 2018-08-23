import json
import os

# This loads the json file and stores the amino acid exact masses as a dictionary
def getAminoAcidMasses():
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    with open(os.path.join(__location__, 'aaIonMasses.json'), encoding='utf8') as f:
        aaIonMasses = json.load(f)
    return aaIonMasses



# This function makes the sequence file
def createSequenceFile():
    seqFile = "sequence.txt"
    f = open(seqFile, 'w')
    f.write(input("\nEnter your sequence using the availible codes, seperated by spaces...\n\n"))
    f.close()
    print("\n")
    return


# This function will take the string from the sequence file and make the list of amino acids
def getSequenceAAs():

    dictOfAA = getAminoAcidMasses()
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    with open(os.path.join(__location__, 'sequence.txt'), encoding='utf8') as f:
        seq = f.readline().rstrip().split(' ')

    cleanedSequence = [x for x in seq if x in dictOfAA.keys()]
    badEntries = [x for x in seq if x not in dictOfAA.keys()]

    # Warns the user if there is a possible typo in the sequence
    if badEntries:
        print('Your seuqnece had amino acids for which I do not know the mass')
        print('You should fix your sequence file')
        print('I removed the following amino acids from the sequence')
        for AA in badEntries:
            print(AA)
        print("")

    return cleanedSequence


# This funciton accepts a list of amino acids and returns the exact mass of that sequence
def getExactMass(listOfAA, dictOfAA):
    exactMass = 0.0
    for i in listOfAA:
        exactMass += dictOfAA[i]
    return exactMass

# This function accepts the sequence list, and returns a list that has the unique amino acids
# i.e. no repetition...
def getUniqueAA(sequence):
    uniqueAA = []
    for i in sequence:
        if i not in uniqueAA:
            uniqueAA.append(i)

    return uniqueAA

# This function builds the dictionary ofcalculated masses of ions,
# using same process as the iOS app
def populateMassDictionary(mass, sequence, uniqueAA, dictOfAA):
    massDictionary = {} # Formating will be {Float : String}
    # exact masses of the ions used throughout script
    H = 1.00782; Na = 22.9898; K = 38.96371; O = 15.9996; dS = 2.01564

    # TODO let the user specify these at some point
    maxCharge = 3
    maxAggregate = 2
    z = 0

    # define all of the variables to keep track of things
    numTFAester = 0
    numBoc = 0
    numtBu = 0
    numTrt = 0
    numGuan = 0

    # counts how many, and what kind, of side reactions can occur
    for i in sequence:
        if i == "lys" or i == "orn" or i == "b-lys":
            numBoc += 1
            numGuan += 1
        elif i == "ser" or i == "thr" or i == "b-ser" or i == "b-thr":
            numtBu += 1
            numTFAester += 1
        elif i =="his" or i == "b-his":
            numTFAester +=1
        elif i == "cys" or i == "b-cys":
            numTrT += 1
        elif i == "HO" or i == "NH2" or i == "NHCH3":
            numGuan += 1
        else:
            continue

    # This is the loop that loops through the possible ions, adds them to the dictionary
    # l is index on M
    for l in range(1, maxAggregate + 1):
        # i is index of H
        for i in range(0, maxCharge + 1):
            # j is index of Na
            for j in range(0, maxCharge + 1):
                # k is index of H
                for k in range(0, maxCharge + 1):
                    ionMass = mass * float(l) + (float(i) * H + float(j) * Na + float(k) * K )

                    if i == 0 and j == 0 and k == 0:
                        continue
                    elif i + j + k > maxCharge:
                        continue

                    # H masses
                    elif i >= 1 and j == 0 and k == 0:
                        ionMass = ionMass / float(i)
                        z = i
                        massDictionary[ionMass] = f"[{l}M + {i}H]{z}+"

                    # Na masses
                    elif i == 0 and j >= 1 and k == 0:
                        ionMass = ionMass / float(j)
                        z = j
                        massDictionary[ionMass] = f"[{l}M + {j}Na]{z}+"

                    # K masses
                    elif i == 0 and j == 0 and k >= 1:
                        ionMass = ionMass / float(k)
                        z = k
                        massDictionary[ionMass] = f"[{l}M + {k}K]{z}+"

                    # H Na masses
                    elif i >= 1 and j >= 0 and k == 0:
                        ionMass = ionMass / (float(i) + float(j))
                        z = i + j
                        massDictionary[ionMass] = f"[{l}M + {i}H + {j}Na]{z}+"

                    # H K masses
                    elif i >= 1 and j == 0 and k >= 0:
                        ionMass = ionMass / (float(i) + float(k))
                        z = i + k
                        massDictionary[ionMass] = f"[{l}M + {i}H + {k}K]{z}+"

                    # Na K masses
                    elif i == 0 and j >= 0 and k >= 0:
                        ionMass = ionMass / (float(j) + float(k))
                        z = j + k
                        massDictionary[ionMass] = f"[{l}M + {j}Na + {k}K]{z}+"

                    # H Na K masses
                    elif i >= 1 and j >= 0 and k >= 0:
                        ionMass = ionMass / (float(i) + float(j) + float(k))
                        z = i + j + k
                        massDictionary[ionMass] = f"[{l}M + {i}H + {j}Na + {k}K]{z}+"
                    else:
                        continue

                    # Calculates Boc additions, if its possible
                    if numBoc > 0:
                        for x in range(1, numBoc + 1):
                            ionMass = (mass + dictOfAA['Boc']) * float(l) + (float(i) * H + float(j) * Na + float(k) * K )

                            if i == 0 and j == 0 and k == 0:
                                continue
                            elif i + j + k > maxCharge:
                                continue

                            # H masses
                            elif i >= 1 and j == 0 and k == 0:
                                ionMass = ionMass / float(i)
                                z = i
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {i}H]{z}+"

                            # Na masses
                            elif i == 0 and j >= 1 and k == 0:
                                ionMass = ionMass / float(j)
                                z = j
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {j}Na]{z}+"

                            # K masses
                            elif i == 0 and j == 0 and k >= 1:
                                ionMass = ionMass / float(k)
                                z = k
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {k}K]{z}+"

                            # H Na masses
                            elif i >= 1 and j >= 0 and k == 0:
                                ionMass = ionMass / (float(i) + float(j))
                                z = i + j
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {i}H + {j}Na]{z}+"

                            # H K masses
                            elif i >= 1 and j == 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(k))
                                z = i + k
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {i}H + {k}K]{z}+"

                            # Na K masses
                            elif i == 0 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(j) + float(k))
                                z = j + k
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {j}Na + {k}K]{z}+"

                            # H Na K masses
                            elif i >= 1 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(j) + float(k))
                                z = i + j + k
                                massDictionary[ionMass] = f"[{l}M + {x}Boc + {i}H + {j}Na + {k}K]{z}+"
                            else:
                                continue

                    # Calculates tBu
                    if numtBu > 0:
                        for x in range(1, numtBu + 1):
                            ionMass = (mass + dictOfAA['tBu']) * float(l) + (float(i) * H + float(j) * Na + float(k) * K )

                            if i == 0 and j == 0 and k == 0:
                                continue
                            elif i + j + k > maxCharge:
                                continue

                            # H masses
                            elif i >= 1 and j == 0 and k == 0:
                                ionMass = ionMass / float(i)
                                z = i
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {i}H]{z}+"

                            # Na masses
                            elif i == 0 and j >= 1 and k == 0:
                                ionMass = ionMass / float(j)
                                z = j
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {j}Na]{z}+"

                            # K masses
                            elif i == 0 and j == 0 and k >= 1:
                                ionMass = ionMass / float(k)
                                z = k
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {k}K]{z}+"

                            # H Na masses
                            elif i >= 1 and j >= 0 and k == 0:
                                ionMass = ionMass / (float(i) + float(j))
                                z = i + j
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {i}H + {j}Na]{z}+"

                            # H K masses
                            elif i >= 1 and j == 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(k))
                                z = i + k
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {i}H + {k}K]{z}+"

                            # Na K masses
                            elif i == 0 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(j) + float(k))
                                z = j + k
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {j}Na + {k}K]{z}+"

                            # H Na K masses
                            elif i >= 1 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(j) + float(k))
                                z = i + j + k
                                massDictionary[ionMass] = f"[{l}M + {x}tBu + {i}H + {j}Na + {k}K]{z}+"
                            else:
                                continue

                    # Calculates Trt and Acm
                    if numTrt > 0:
                        for x in range(1, numTrT + 1):
                            ionMass = (mass + dictOfAA['Trt']) * float(l) + (float(i) * H + float(j) * Na + float(k) * K )
                            ionMass1 = (mass + dictOfAA['Acm']) * float(l) + (float(i) * H + float(j) * Na + float(k) * K )

                            if i == 0 and j == 0 and k == 0:
                                continue
                            elif i + j + k > maxCharge:
                                continue

                            # H masses
                            elif i >= 1 and j == 0 and k == 0:
                                ionMass = ionMass / float(i)
                                z = i
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {i}H]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {i}H]{z}+"

                            # Na masses
                            elif i == 0 and j >= 1 and k == 0:
                                ionMass = ionMass / float(j)
                                z = j
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {j}Na]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {j}Na]{z}+"

                            # K masses
                            elif i == 0 and j == 0 and k >= 1:
                                ionMass = ionMass / float(k)
                                z = k
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {k}K]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {k}K]{z}+"

                            # H Na masses
                            elif i >= 1 and j >= 0 and k == 0:
                                ionMass = ionMass / (float(i) + float(j))
                                z = i + j
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {i}H + {j}Na]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {i}H + {j}Na]{z}+"

                            # H K masses
                            elif i >= 1 and j == 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(k))
                                z = i + k
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {i}H + {k}K]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {i}H + {k}K]{z}+"

                            # Na K masses
                            elif i == 0 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(j) + float(k))
                                z = j + k
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {j}Na + {k}K]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {j}Na + {k}K]{z}+"

                            # H Na K masses
                            elif i >= 1 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(j) + float(k))
                                z = i + j + k
                                massDictionary[ionMass] = f"[{l}M + {x}Trt + {i}H + {j}Na + {k}K]{z}+"
                                massDictionary[ionMass1] = f"[{l}M + {x}Acm + {i}H + {j}Na + {k}K]{z}+"
                            else:
                                continue

                    # Calculates TFA esters
                    if numTFAester > 0:
                        for x in range(1, numTFAester + 1):
                            ionMass = (mass + dictOfAA['TFAester']) * float(l) + (float(i) * H + float(j) * Na + float(k) * K )

                            if i == 0 and j == 0 and k == 0:
                                continue
                            elif i + j + k > maxCharge:
                                continue

                            # H masses
                            elif i >= 1 and j == 0 and k == 0:
                                ionMass = ionMass / float(i)
                                z = i
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {i}H]{z}+"

                            # Na masses
                            elif i == 0 and j >= 1 and k == 0:
                                ionMass = ionMass / float(j)
                                z = j
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {j}Na]{z}+"

                            # K masses
                            elif i == 0 and j == 0 and k >= 1:
                                ionMass = ionMass / float(k)
                                z = k
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {k}K]{z}+"

                            # H Na masses
                            elif i >= 1 and j >= 0 and k == 0:
                                ionMass = ionMass / (float(i) + float(j))
                                z = i + j
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {i}H + {j}Na]{z}+"

                            # H K masses
                            elif i >= 1 and j == 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(k))
                                z = i + k
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {i}H + {k}K]{z}+"

                            # Na K masses
                            elif i == 0 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(j) + float(k))
                                z = j + k
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {j}Na + {k}K]{z}+"

                            # H Na K masses
                            elif i >= 1 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(j) + float(k))
                                z = i + j + k
                                massDictionary[ionMass] = f"[{l}M + {x}TFA esters + {i}H + {j}Na + {k}K]{z}+"
                            else:
                                continue

                    # Calcualtes gunidinylations
                    if numGuan > 0:
                        for x in range(1, numGuan + 1):
                            ionMass = (mass + dictOfAA['Guan']) * float(l) + (float(i) * H + float(j) * Na + float(k) * K )

                            if i == 0 and j == 0 and k == 0:
                                continue
                            elif i + j + k > maxCharge:
                                continue

                            # H masses
                            elif i >= 1 and j == 0 and k == 0:
                                ionMass = ionMass / float(i)
                                z = i
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {i}H]{z}+"

                            # Na masses
                            elif i == 0 and j >= 1 and k == 0:
                                ionMass = ionMass / float(j)
                                z = j
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {j}Na]{z}+"

                            # K masses
                            elif i == 0 and j == 0 and k >= 1:
                                ionMass = ionMass / float(k)
                                z = k
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {k}K]{z}+"

                            # H Na masses
                            elif i >= 1 and j >= 0 and k == 0:
                                ionMass = ionMass / (float(i) + float(j))
                                z = i + j
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {i}H + {j}Na]{z}+"

                            # H K masses
                            elif i >= 1 and j == 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(k))
                                z = i + k
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {i}H + {k}K]{z}+"

                            # Na K masses
                            elif i == 0 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(j) + float(k))
                                z = j + k
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {j}Na + {k}K]{z}+"

                            # H Na K masses
                            elif i >= 1 and j >= 0 and k >= 0:
                                ionMass = ionMass / (float(i) + float(j) + float(k))
                                z = i + j + k
                                massDictionary[ionMass] = f"[{l}M + {x}Guanidinylation + {i}H + {j}Na + {k}K]{z}+"
                            else:
                                continue
    return massDictionary

# This function takes in the queryed mass, and the mass dictionary,
# searches for a match within that dictionary
# returns a list of matches
def getMassMatch(query, possibleMasses):
    matches = []
    for mass, ion in possibleMasses.items():
        if abs(query - mass) <= 0.3:
            matches.append(ion + ' : ' + str(round(mass, 4)))
        else:
            continue

    if not matches:
        matches.append("No matches found")

    return matches


# This function writes all possible masses to a text file
def writeMassDictionary(possibleMasses):
    massFile = "possible masses.txt"
    f = open(massFile, 'w')
    for mass, ion in possibleMasses.items():
        string = ion +' : ' + str(round(mass, 4))
        f.write(string)
        f.write('\n')
    f.close()

# For testing purposes
if __name__ == '__main__':
    quit()
