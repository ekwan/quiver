# This file parses output from QUIVER and calculates the KIEs.
# Eugene Kwan, November 2015
#
# usage: awk -f quiver_analysis.awk gs.q3 gs.qout ts.qout

BEGIN {
    print "\n=== QUIVER Analysis ===\n"
}

# reset variables at the beginning of the file
FNR == 1 {
    fileCount++
    partitionCount = 0
    frequencyCount = 0
    filenames[fileCount]=FILENAME
}

# read the temperature
fileCount == 1 && FNR == 1 {
    temperature = $1
    printf "Temperature: %.f K\n", temperature
}

# KIEs will be referenced to this isotopologue
fileCount == 1 && FNR == 2 {
    referenceIsotopologue = $1
    if ( referenceIsotopologue == 0 )
        print "KIEs are not referenced to any nucleus."
    else
        printf "KIEs are referenced to isotopologue number %d.\n", referenceIsotopologue
}

# read the description of each isotopologue
fileCount == 1 && FNR > 2 {
    descriptions[FNR-2]=$0
    numberOfIsotopologues=FNR-2
}

# read the reduced isotopic partition functions
/S2/ {
    partitionCount++
    partitionFunctions[fileCount,partitionCount]=$2
    # print fileCount, partitionCount, $2
}

# read the transition state imaginary frequencies
/FREQ/ && /    -/ && fileCount == 3 {
    frequencyCount++
    frequencies[frequencyCount]=$3
    #print fileCount, frequencyCount, $3
}

# calculate the KIEs
END {
    print "\nNote: isotopomer descriptions refer to ground state atom numbers."

    # calculate some quantities we need for the tunnelling corrections
    if ( frequencyCount > 0 )
        {
            print "\nThis is a kinetic isotope effect calculation.  Tunnelling corrections have been applied,"
            print "but may not be accurate for H/D/T KIEs."
            u_unsubstituted = -1.43877 * frequencies[1] / temperature
            u_unsubstituted_temp = u_unsubstituted / sin(u_unsubstituted*0.5)
            temp1 = 1 + u_unsubstituted*u_unsubstituted/24
            for (i=1; i <= numberOfIsotopologues; i++)
                {
                    S2overS1 = partitionFunctions[2,i]/partitionFunctions[3,i]

                    # this is the KIE with no tunnelling
                    rawKIE = S2overS1 * frequencies[1] / frequencies[i+1]

                    # this is the Bell infinite parabola correction
                    u_substituted = -1.43877 * frequencies[i+1] / temperature
                    u_substituted_temp = sin(u_substituted*0.5) / u_substituted
                    bellCorrection = u_unsubstituted_temp * u_substituted_temp
                    infiniteParabolaKIE = rawKIE * bellCorrection
                    
                    # this is the Widmer correction
                    temp2 = 1 + u_substituted*u_substituted/24
                    widmerCorrection = temp1/temp2
                    widmerKIE = rawKIE * widmerCorrection
                    
                    # store uncorrected KIEs
                    KIE[i,1]=rawKIE
                    KIE[i,2]=infiniteParabolaKIE
                    KIE[i,3]=widmerKIE
                }
        }
    else
        {   
            print "\nThis is an equilibrium isotope effect calculation.  No tunnelling corrections have been applied."
            for (i=1; i <= numberOfIsotopologues; i++)
                {
                    rawKIE = partitionFunctions[2,i]/partitionFunctions[3,i]
            
                    # store uncorrected KIEs
                    KIE[i,1]=rawKIE
                    KIE[i,2]=rawKIE
                    KIE[i,3]=rawKIE
                }
        }
    
    printf "\n"
    print "isotopologue description                                  uncorrected      Widmer     infinite parabola"
    print "                                                              KIE           KIE              KIE"
    
    referenceIsotopologue == 0 ? referenceKIE[1] = 1.000 : referenceKIE[1] = KIE[referenceIsotopologue,1]
    referenceIsotopologue == 0 ? referenceKIE[2] = 1.000 : referenceKIE[2] = KIE[referenceIsotopologue,2]
    referenceIsotopologue == 0 ? referenceKIE[3] = 1.000 : referenceKIE[3] = KIE[referenceIsotopologue,3]
    
    if ( referenceIsotopologue > 0 )
        printf "%60-s %5.3f         %5.3f            %5.3f\n", "reference KIEs", referenceKIE[1], referenceKIE[3], referenceKIE[2]
    for (i=1; i <= numberOfIsotopologues; i++)
        {
            # print result
            rawKIE = KIE[i,1] / referenceKIE[1]
            widmerKIE = KIE[i,3] / referenceKIE[3]
            infiniteParabolaKIE = KIE[i,2] / referenceKIE[2]
            if ( i != referenceIsotopologue )
                printf "%60-s %5.3f         %5.3f            %5.3f\n", descriptions[i], rawKIE, widmerKIE, infiniteParabolaKIE
            else
                printf "%60-s %5.3f         %5.3f            %5.3f\n", "[" descriptions[i] "] (referenced to 1.00)", rawKIE, widmerKIE, infiniteParabolaKIE
        }
}
