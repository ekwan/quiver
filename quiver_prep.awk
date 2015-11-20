# This program prepares input files for QUIVER.
# Eugene Kwan, November 2015
#
# usage: awk -f quiver_prep.awk g09_gs_freq.out g09_ts_freq.out
#
# - frequency files should have #p freq specified
# - control all other program behavior in the BEGIN section below

BEGIN {

# configure the behavior of QUIVER here
scaling = 0.9614                    # scaling factor for frequencies, should be adjusted for whatever level of theory is being used
temperature = 273+120               # in K

# this program will ask QUIVER to calculate the KIEs the single-atom isotope replacements given below
# for example, replacements[6]=13 would replace 12 C1/12 C2/12 C2 with three isotopomers,
# 13 C1 / 12 C2 / 12 C3
# 12 C1 / 13 C2 / 12 C3
# 12 C1 / 12 C2 / 13 C3
replacements[6]=13                  # replace carbon atoms (atomic number 6) with carbon-13s
replacements[8]=17                  # replace oxygen atoms (atomic number 8) with oxygen-17s

# calculate a special isotopologue if desired
special = 0                         # set to 1 if you want to create a custom isotopologue
specialIsotopes[7]=2                # as an example, this would replace atoms 7 and 8 with deuteriums
specialIsotopes[8]=2

# KIEs will be referenced to this isotopologue
# 0 : do not reference KIEs
# 1 : assume the KIE in the first isotopologue is 1.000, and divide all others by the predicted KIE there
# ...and so forth
# To figure out which isotopologue to use, run this first with referenceIsotopologue set to zero,
# then look at the number of the isotopologue in the output.
referenceIsotopologue = 5

# store some data about atomic numbers and their standard weights
atomicWeight[1]=1     ; atomicSymbol[1] ="H"     # hydrogen
atomicWeight[6]=12    ; atomicSymbol[6] ="C"     # carbon
atomicWeight[7]=14    ; atomicSymbol[7] ="N"     # nitrogen
atomicWeight[8]=16    ; atomicSymbol[8] ="O"     # oxygen
atomicWeight[9]=19    ; atomicSymbol[9] ="F"     # fluorine
atomicWeight[16]=32   ; atomicSymbol[16]="S"     # sulfur
atomicWeight[17]=35   ; atomicSymbol[17]="Cl"    # chlorine
atomicWeight[22]=48   ; atomicSymbol[22]="Ti"    # titanium
atomicWeight[35]=80   ; atomicSymbol[35]="Br"    # bromine
}

########## Read Data Files #########

# this is the start of any file
FNR == 1 {
    fileCount++
    filenames[fileCount]=FILENAME
}

# read the title of each file
/l101.exe/,/Symbolic Z-matrix/ {
    if ( match($0,"------") > 0 && length(title[fileCount]) == 0 )
        {
            getline
            while ( match($0,"------") == 0 )
                {
                    title[fileCount] = title[fileCount] trim($0)
                    getline
                }
        }
}

# read the last geometry of each file
/Standard orientation/,/Rotational constants/ {
    if ( NF == 6 && match($0,"[a-zA-Z]") == 0 )
        {
            x_coordinates[fileCount,$1] = $4
            y_coordinates[fileCount,$1] = $5
            z_coordinates[fileCount,$1] = $6
            atomicNumbers[fileCount,$1] = $2
            numberOfAtoms[fileCount] = $1
        }
}

# read the number of imaginary frequencies
/imaginary frequencies ignored/ {
    numberOfImaginaries[fileCount] = $1
}

# read the cartesian force constants from the archive
/l9999.exe/,/\\\@/ {
    archive[fileCount] = archive[fileCount] trim($0)
}

########## Print Input Deck #########

# check the inputs carefully and calculate the KIEs
END {

# put in default titles if necessary
if ( length(title[1]) == 0 )
    title[1] = "ground state"
if ( length(title[2]) == 0 )
    title[2] = "transition state"
printf "Read ground state file from %s (%s).\n", filenames[1], title[1]
printf "Read transition state file from %s (%s).\n", filenames[2], title[2]

# check geometries are present and consistent
if ( numberOfAtoms[1] != numberOfAtoms[2] )
    abort(sprintf("First structure has %d atoms, but second structure has %d atoms!\n", numberOfAtoms[1], numberOfAtoms[2]))
if ( numberOfAtoms[1] == 0 || numberOfAtoms[2] == 0 ) abort("No geometry read for one of the structures.")
for (i=1; i <= numberOfAtoms[1]; i++)
    {
        if ( atomicNumbers[1,i] != atomicNumbers[2,i] )
        {
            print "Mismatch in atomic numbers between files!  Files must have the same atom ordering."
            printf "Discrepancy for atom %d: atomic numbers are %d and %d\n", i, atomicNumbers[1,i], atomicNumbers[2,i]
            exit 1
        }
    }
printf "%d atoms read successfully.\n", numberOfAtoms[1]

# check imaginary frequencies
printf "%d imaginary frequencies in ground state file.\n", numberOfImaginaries[1]
printf "%d imaginary frequencies in transition state file.\n", numberOfImaginaries[2]

# parse the cartesian force constants
cartesianForceConstants[1] = getForceConstants(archive[1])
cartesianForceConstants[2] = getForceConstants(archive[2])
printf "%d cartesian force constants read from ground state file.\n", newlines(cartesianForceConstants[1]) 
printf "%d cartesian force constants read from transition state file.\n", newlines(cartesianForceConstants[2])

# print the output files
printQuiverInput(1, "gs")
printQuiverInput(2, "ts")

# tell the user we are done
print "QUIVER input decks placed in gs.q1 and ts.q1."
print "Cartesian force constants placed in gs.q1 and ts.q1."
print "QUIVER will need the q1 and q2 files renamed to temp.q1 and temp.q2 to run."
print "Descriptions of all isotopomers written to gs.q3 and ts.q3, which are not needed by QUIVER."
print "\n*** QUIVER files written successfully! ***"
}

########## Helper Functions #########

# quit with the given message
function abort(string) {
    print string
    exit 1
}

# functions to trim whitespace
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }

# count number of newlines
function newlines(s) {
    return split(s, tempfields, "\n") - 1
}

########## Parser Functions #########

# extracts the cartesian force constants from a string
function getForceConstants(archiveString) {
    numberOfFields = split(archiveString, fields, "\\")
    
    targetFieldNumber = -1
    for (i=1; i <= numberOfFields; i++)
        {
            if ( match(tolower(fields[i]), "nimag=") > 0 )
                {
                    targetFieldNumber = i + 2
                    break
                }
        }
    if ( targetFieldNumber == -1 ) abort("Problem reading Cartesian force constants from archive!")
    
    forceConstantString = ""
    numberOfFields = split(fields[targetFieldNumber], fields, ",")
    if ( numberOfFields == 0 ) abort("Problem reading Cartesian force constants from archive!")
    for (i=1; i <= numberOfFields; i++)
        forceConstantString = forceConstantString trim(fields[i]) "\n"
    return forceConstantString
}

# prints the QUIVER q1 and q2 input files
function printQuiverInput(fileIndex, prefix) {
    # print the main input deck 
    outputFilename = prefix ".q1"
    printf "" > outputFilename
    printf "title: %s\n", title[fileIndex] >> outputFilename

    # this is the number of isotopologues that will be computed
    # note that this includes the unsubstituted isotopologue
    numberOfReplacements = 0
    for (i=1; i <= numberOfAtoms[fileIndex]; i++)
        if (atomicNumbers[fileIndex,i] in replacements)
            numberOfReplacements++
    if ( ! (special == 0 || special == 1) ) abort("check value of special in BEGIN block")
    numberOfIsotopologues = 1 + numberOfReplacements + special
    print "  " numberOfIsotopologues >> outputFilename

    # print the geometry
    print "parent" >> outputFilename
    print "  " numberOfAtoms[fileIndex] >> outputFilename
    for (i=1; i <= numberOfAtoms[fileIndex]; i++)
        printf "   %12.6f%12.6f%12.6f\n", x_coordinates[fileIndex,i], y_coordinates[fileIndex,i], z_coordinates[fileIndex,i] >> outputFilename

    # print the standard isotopologue
    for (i=1; i <= numberOfAtoms[fileIndex]; i++)
        {
            thisAtomicNumber = atomicNumbers[fileIndex,i]
            if ( ! thisAtomicNumber in atomicWeight )
                abort(sprintf("unrecognized atomic number %d!\n", thisAtomicNumber))
            printf atomicWeight[thisAtomicNumber] >> outputFilename
            if ( i < numberOfAtoms[fileIndex] )
                printf "," >> outputFilename
        }
    printf "\n" >> outputFilename

    # print the scaling factor and temperature
    print scaling >> outputFilename
    print "   1" >> outputFilename
    print temperature >> outputFilename

    # print the isotopologues for the standard replacements (e.g., 12C with 13C, one at a time)
    replacementCount = 0
    for (i=1; i <= numberOfAtoms[fileIndex]; i++)
        {
            replaced = 0
            isotopologueString1 = ""
            isotopologueString2 = ""
            for (j=1; j <= numberOfAtoms[fileIndex]; j++)
                {
                    thisAtomicNumber = atomicNumbers[fileIndex,j]
                    if ( i==j && thisAtomicNumber in replacements )
                        {
                            # this is a replacement
                            isotopologueString1 = sprintf("%d%s @ atom %d", replacements[thisAtomicNumber],atomicSymbol[thisAtomicNumber],j)
                            isotopologueString2 = isotopologueString2 replacements[thisAtomicNumber]
                            replaced = 1
                        }
                    else
                        isotopologueString2 = isotopologueString2 atomicWeight[thisAtomicNumber]
                    if ( j < numberOfAtoms[fileIndex] )
                        isotopologueString2 = isotopologueString2 ","
                }
            if ( replaced == 1 )
                {
                    replacementCount++
                    # this is the name of the isotopomer, can be anything
                    descriptions[replacementCount] = isotopologueString1
                    printf "%d - %s\n", replacementCount, isotopologueString1 >> outputFilename
                    print isotopologueString2 >> outputFilename
                }
        }

    # print the special isotopologue if requested
    if ( special == 1 ) 
        {
            # check the input is sensible
            for ( i in specialIsotopes )
                if ( i+0 < 1 || i+0 > numberOfAtoms[fileIndex]+0 )
                    abort(sprintf("special isotope atom number of of range: %d (allowed range 1-%d)", i, numberOfAtoms[fileIndex]))
            # make the special isotopologue
            isotopologueString1 = "special"
            isotopologueString2 = ""
            for (i=1; i <= numberOfAtoms[fileIndex]; i++)
                {
                    thisAtomicNumber = atomicNumbers[fileIndex,i]
                    if ( i in specialIsotopes )
                        isotopologueString2 = isotopologueString2 specialIsotopes[i]
                    else
                        isotopologueString2 = isotopologueString2 atomicWeight[thisAtomicNumber]
                    if ( i < numberOfAtoms[fileIndex] )
                        isotopologueString2 = isotopologueString2 ","
                }
            print isotopologueString1 >> outputFilename
            print isotopologueString2 >> outputFilename
        }

    # print the summary of changes
    printf " %d\n", numberOfReplacements+special >> outputFilename
    for (i=1; i <= numberOfReplacements+special; i++)
        printf "   1  %d\n", i+1 >> outputFilename

    # print the cartesian force constants
    outputFilename = prefix ".q2"
    print cartesianForceConstants[fileIndex] > outputFilename

    # check reference isotopologue
    if ( referenceIsotopologue < 0 || referenceIsotopologue > numberOfReplacements + special )
        abort("check reference isotopologue")

    # print a list of all the descriptions of what each isotopologue is
    outputFilename = prefix ".q3"
    print temperature > outputFilename
    print referenceIsotopologue >> outputFilename
    for (i=1; i <= replacementCount; i++)
        print descriptions[i] >> outputFilename
    if ( special == 1 )
        print "special" >> outputFilename
}
