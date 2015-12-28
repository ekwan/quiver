# This program prepares input files for QUIVER.
# Eugene Kwan, November 2015
#
# usage: awk -f quiver_prep.awk quiver.config g09_gs_freq.out g09_ts_freq.out
#
# - frequency files should have #p freq specified

BEGIN {
# store some data about atomic numbers and their standard weights
# these weights will be used in the "unsubstituted" isotopomer
atomicWeight[1]=1     ; atomicSymbol[1] ="H"     # hydrogen
atomicWeight[6]=12    ; atomicSymbol[6] ="C"     # carbon
atomicWeight[7]=14    ; atomicSymbol[7] ="N"     # nitrogen
atomicWeight[8]=16    ; atomicSymbol[8] ="O"     # oxygen
atomicWeight[9]=19    ; atomicSymbol[9] ="F"     # fluorine
atomicWeight[16]=32   ; atomicSymbol[16]="S"     # sulfur
atomicWeight[17]=35   ; atomicSymbol[17]="Cl"    # chlorine
atomicWeight[22]=48   ; atomicSymbol[22]="Ti"    # titanium
atomicWeight[35]=80   ; atomicSymbol[35]="Br"    # bromine

# default values, will be overriden by config file
scaling = 1.00
temperature = 273
referenceIsotopomer = 0
}

########## Read Data Files #########

# this is the start of any file
NR != FNR && FNR == 1 {
    fileCount++
    filenames[fileCount]=FILENAME
}

# read the scaling factor
NR == FNR {
    fieldName = tolower($1)
    if ( fieldName == "scaling" )
        {
            if ( NF != 2 )
                abort("check config line (unexpected number of fields)\n" $0)
            if ( $2 < 0.8 || $2 > 1.2 )
                abort("check config line (unusual scaling factor)\n" $0)
            scaling = $2
        }
    else if ( fieldName == "temperature" )
        {
            if ( NF != 2 )
                abort("check config line (unexpected number of fields)\n" $0)
            if ( $2 < 0 )
                abort("check config line (negative temperature)\n" $0)
            temperature = $2
        }
    else if ( fieldName == "isotopomer" )
        {
            if ( NF != 5 )
                abort("check config line (unexpected number of fields)\n" $0)
            if ( $2 < 1 || $3 < 1 || $4 < 1 || $5 < 1 )
                abort("check config line (unexpected atom number):\n" $0)
            if ( $2 < numberOfIsotopomers || $2 > numberOfIsotopomers + 1 )
                abort("check config line (non-sequential isotopomer number):\n" $0)
            if ( numberOfIsotopomers == 0 && $2 != 1 )
                abort("check config line (must start at isotopomer 1):\n" $0)
            numberOfIsotopomers=$2
            isotopomerLines++
            isotopomerEntries[isotopomerLines]=$2 " " $3 " " $4 " " $5
        }
    else if ( fieldName == "reference_isotopomer" )
        {
            if ( NF != 2 )
                abort("check config line (unexpected number of fields)\n" $0)
            if ( $2 < 0 )
                abort("negative reference isotopomer not allowed")
            referenceIsotopomer = $2
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

# check atomic weights and symbols are available
for (i=1; i <= 2; i++)
    {
        for (j=1; j <= numberOfAtoms[i]; j++)
            {
                thisAtomicNumber = atomicNumbers[i,j]
                if (! thisAtomicNumber in atomicWeight)
                    abort("no known atomic weight for Z=" thisAtomicNumber "!")
                if (! thisAtomicNumber in atomicSymbol)
                    abort("no known atomic symbol for Z=" thisAtomicNumber "!")
            }
    }

# put in default titles if necessary
title[1] = "ground state " filenames[1]
title[2] = "transition state " filenames[2]
printf "Read ground state file from %s (%d atoms).\n", filenames[1], numberOfAtoms[1]
printf "Read transition state file from %s (%d atoms).\n", filenames[2], numberOfAtoms[2]

# check imaginary frequencies
printf "%d imaginary frequencies in ground state file.\n", numberOfImaginaries[1]
printf "%d imaginary frequencies in transition state file.\n", numberOfImaginaries[2]

# parse the cartesian force constants
print "Reading force constants..."
cartesianForceConstants[1] = getForceConstants(archive[1])
cartesianForceConstants[2] = getForceConstants(archive[2])
printf "%d cartesian force constants read from ground state file.\n", newlines(cartesianForceConstants[1]) 
printf "%d cartesian force constants read from transition state file.\n", newlines(cartesianForceConstants[2])

# check isotopomers are sensible
if (numberOfIsotopomers < 1)
    abort("Must specify at least one isotopomer!")
for (i=1; i <= isotopomerLines; i++)
    {
        split(isotopomerEntries[i],fields," ")
        thisLine = "isotopomer " fields[1] " " fields[2] " " fields[3] " " fields[4]
        
        gsAtomNumber = fields[2]
        tsAtomNumber = fields[3]
        replacementWeight = fields[4]

        # check the gs and ts atom numbers are in range for their respective geometries
        if ( gsAtomNumber > numberOfAtoms[1] || gsAtomNumber < 1 )
            abort("ground state atom number " gsAtomNumber " out of range for\n" thisLine)
        if ( tsAtomNumber > numberOfAtoms[2] || tsAtomNumber < 1 )
            abort("transition state atom number " tsAtomNumber " out of range for\n" thisLine)

        # check that the gs and ts atom numbers correspond to the same atomic number
        # that is, they might have different atom numbers, but they should both be replacing the same kind of atom
        if ( atomicNumbers[1,gsAtomNumber] != atomicNumbers[2,tsAtomNumber] )
            abort("atomic numbers do not match for this line:\n" thisLine)

        # check that the replacement weight is reasonable
        thisAtomicNumber = atomicNumbers[1,gsAtomNumber]
        normalWeight = atomicWeight[thisAtomicNumber]
        if ( replacementWeight > normalWeight + 5 || replacementWeight < normalWeight - 5 )
            abort(sprintf("super light or heavy replacement %s --> %s on line:\n%s", normalWeight, replacementWeight, thisLine))
    }
printf "%d isotopomers read.\n", numberOfIsotopomers

# print the output files
printQuiverInput(1, "gs")
printQuiverInput(2, "ts")

# tell the user we are done
print "QUIVER input decks placed in gs.q1 and ts.q1."
print "Cartesian force constants placed in gs.q2 and ts.q2."
print "Summary data placed in gs.q3 and ts.q3."
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

    # write the number of isotopomers that will be computed by quiver
    # note that this includes the unsubstituted isotopomer
    print "  " (numberOfIsotopomers+1) >> outputFilename

    # print the number of atoms
    print "parent" >> outputFilename
    print "  " numberOfAtoms[fileIndex] >> outputFilename
    
    # print the geometry
    for (i=1; i <= numberOfAtoms[fileIndex]; i++)
        printf "   %12.6f%12.6f%12.6f\n", x_coordinates[fileIndex,i], y_coordinates[fileIndex,i], z_coordinates[fileIndex,i] >> outputFilename

    # print the standard isotopomer
    for (i=1; i <= numberOfAtoms[fileIndex]; i++)
        {
            thisAtomicNumber = atomicNumbers[fileIndex,i]
            printf atomicWeight[thisAtomicNumber] >> outputFilename
            if ( i < numberOfAtoms[fileIndex] )
                printf "," >> outputFilename
        }
    printf "\n" >> outputFilename

    # print the scaling factor and temperature
    print scaling >> outputFilename
    print "   1" >> outputFilename       # number of temperatures is always 1
    print temperature >> outputFilename

    # loop over each replacement isotopomer
    delete descriptions
    for (i=1; i <= numberOfIsotopomers; i++)
        {
            # grab all replacements for this isotopomer
            delete replacement
            for (j=1; j <= isotopomerLines; j++)
                {
                    split(isotopomerEntries[j],fields," ")
                    if ( fields[1] == i )
                        {
                            if ( fileIndex == 1 )
                                fromAtomNumber = fields[2]
                            else if ( fileIndex == 2 )
                                fromAtomNumber = fields[3]
                            else
                                abort("impossible")

                            if ( fromAtomNumber in replacement )
                                abort("duplicate replacement for isotopomer\n" isotopomerEntries[j])
                            targetWeight = fields[4]
                            replacement[fromAtomNumber] = targetWeight
                        }
                }

            # create a description
            description = "isotopomer " i " -"
            for (j in replacement)
                {
                    thisAtomicNumber = atomicNumbers[fileIndex,j]
                    thisAtomicSymbol = atomicSymbol[thisAtomicNumber]
                    description = description sprintf(" %d%s @ atom %d,", replacement[j], thisAtomicSymbol, j)
                }
            description = substr(description, 1, length(description)-1)

            # write out description
            print description >> outputFilename
            descriptions[i] = description

            # enumerate the atomic weights
            for (j=1; j <= numberOfAtoms[fileIndex]; j++)
                {
                    thisAtomicNumber = atomicNumbers[fileIndex,j]
                    regularWeight = atomicWeight[thisAtomicNumber]
                    if ( j in replacement)
                        printf replacement[j] >> outputFilename
                    else
                        printf regularWeight >> outputFilename
                    if ( j < numberOfAtoms[fileIndex] )
                        printf "," >> outputFilename 
                }
            printf "\n" >> outputFilename
        }

    # print the summary of changes
    printf " %d\n", numberOfIsotopomers >> outputFilename
    for (i=1; i <= numberOfIsotopomers; i++)
        printf "   1  %d\n", i+1 >> outputFilename

    # print the cartesian force constants
    outputFilename = prefix ".q2"
    print cartesianForceConstants[fileIndex] > outputFilename

    # ensure the reference isotopomer is in range
    if ( referenceIsotopomer < 0 || referenceIsotopomer > numberOfIsotopomers )
        abort("check reference isotopologue")

    # print a list of all the descriptions of what each isotopologue is
    outputFilename = prefix ".q3"
    print temperature > outputFilename
    print referenceIsotopomer >> outputFilename
    for (i=1; i <= numberOfIsotopomers; i++)
        print descriptions[i] >> outputFilename
}
