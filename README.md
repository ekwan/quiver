## Scripts for Running QUIVER to Calculate KIEs

This is a series of scripts that facilitate the calculation of kinetic isotope effects using
Bigeleisen-Mayer theory.  You provide Gaussian output files containing the frequencies of the
ground and transition states, as well as which isotopologues you want to compute.  This script
gives the KIEs, with and without standard tunnelling corrections.

As a test case, I have provided the transition states for aliphatic Claisen rearrangement of
allyl vinyl ether that was studied by Singleton in JACS 1999, 121, 10865.  These files are
stored in the `claisen_demo/` folder.  There are two files, one for the ground state (gs) and
one for the transition state (ts).  Both are calculated at B3LYP/6-31G(d).  I lifted the
structures from the SI.  Conveniently, the carbons and oxygen atom are exactly as labeled in
Table 4 of the paper.  The predicted values in the literature for single atom <sup>12</sup>C/<sup>13</sup>C and
<sup>16</sup>O/<sup>17</sup>O substitutions at 120 C are:

```
C1 1.012
C2 0.999
O3 1.017
C4 1.029
C5 1.000 (relative)
C6 1.014
```

Note that this means that the predicted KIEs at all positions are divided by the predicted
KIE at C5.  When I run my script, I get:

```
Temperature: 393 K
KIEs are referenced to isotopologue number 5.

isotopologue description    uncorrected      Widmer     infinite parabola
                                KIE           KIE              KIE
reference KIEs                 1.002         1.002            1.002
13C @ atom 1                   1.011         1.012            1.013
13C @ atom 2                   1.000         1.000            1.000
17O @ atom 3                   1.017         1.018            1.019
13C @ atom 4                   1.028         1.031            1.031
13C @ atom 5                   1.000         1.000            1.000
13C @ atom 6                   1.013         1.015            1.015
```

As you can see, the numbers match to within 0.001.  I'm not exactly sure which tunnelling
correction was used in Table 4.  Singleton notes that tunnelling improves KIE predictions
a bit.  For H/D predictions, though, these crude tunnelling corrections will be poor.  See
the Singleton Claisen paper for a fuller discussion, which references work by Truhlar in
this area.

Below, I explain how to use the scripts.  Note that this uses a scaling factor of 0.9614
for the frequencies.  Consult the references section below for more scaling factors.

# Installation

The first step is to compile QUIVER.  QUIVER is written in Fortran 77 (sigh).  If you are
predicting some unusual elements, it might be necessary to go into the source code and add
the appropriate atomic weight.  The source code is in src/quiver.f.  Use the GNU fortran
compiler (gfortran) to compile it to quiver.exe:

`gfortran quiver.f -o quiver.exe`

If there are any errors, I suggest you check that the whitespace in quiver.f is being
appropriately interpreted on your system.  Copy your compiled version of quiver.exe to
the directory containing the scripts.

# Gaussian Instructions

Using the `#p` directive to request verbose output, obtain frequencies for the ground
and transition state of interest using the `freq` keyword.  (To calculate an equilibrium
isotope effect, replace the transition state with the second ground state.)  The
temperature you choose does not matter and will not be parsed from the file.

I have only tested this in the case where the starting material and transition states
have the same atom numbering.  However, it should work if they are not the same.

# Script Instructions

Run the `quiver.sh` script:

`./quiver.sh claisen_gs.out claisen_ts.out`

The script will run the AWK script quiver_prep.awk to make the QUIVER input files,
run QUIVER, and then run the AWK script quiver_analysis.awk to calculate the KIEs.
The script wil abort if there is a problem (at least in theory).

# Changing Program Behavior

The following discussion references line numbers in `quiver_prep.awk`.

To change the scaling factor for the frequencies, alter line 12.

To change the temperature, alter the temperature on line 13.

To change the reference isotopologue, add the appropriate atom number on line 34.  To
show raw KIEs only, set `referenceIsotopologue` to `0`.

By default, the program will replace any carbon-12 with carbon-13 and any oxygen-16
with oxygen-17, one at a time.  Each isotopologue is then given a number 1, 2, 3, ...
This facilitates the calculation of KIEs for natural
abundance KIE experiments.  See line 15 for an explanation.  To request any other
standard replacements of this type, alter line 20.  It might be necessary to add atomic
weight data to line 37 for any unusual elements.

The program can also calculate a
"special isotopologue."  For example, if you wanted to replace several hydrogens with
deuteriums, you can do so on line 24.  Set special to 1 to calculate a special isotopologue
and indicate the replacements using the `specialIsotopes` array.  In the case of the Claisen,
I have included an example where hydrogens 7 and 8 have been replaced with deuterium.  I
predict an inverse KIE of 0.953 (uncorrected, but referenced to isotopologue 5).

# Technical Details

The scripts extract the last geometry, frequencies, Cartesian force constants, and some
other information from the two g09 .out files for the ground and transition state.  The
main Quiver input is dumped into `gs.q1` and `ts.q1`.  The Cartesian force constants are
placed in `gs.q2` and `ts.q2`.  Descriptions of the isotopologues being calculated are
placed in `gs.q3` and `ts.q3`.  QUIVER does not use the q3 files.  The scripts I wrote use
`gs.q3` (and ignore `ts.q3`) to calculate the KIEs.

As much as possible, I designed the scripts to abort if any part of the parsing didn't
work or got data that make no sense.

# Compatability

This works on my Mac and has only been tested there.  (It is not computataionally demanding
to run--the expensive part of this is the frequencies, which you have to get from Gaussian).
It should work on any Linux system or Cygwin.  I have not tested it with the Windows version
of Gaussian.  Singleton mentions there can be a problem parsing the Cartesian force constants
from the archive at the end of every g09 .out file.  However, I have improved the robustness
of the parsing, so it should always work for regular Linux versions of g09.

The scripts depend on AWK, which should be available just about anywhere.  On Cygwin, be
sure to install the GNU implementation of awk, gawk.

# References

I lifted these references from the Houk group writeup for QUIVER.

Bigeleisen-Mayer theory: J. Chem. Phys.  1947, 15, 261.; Wolfsberg, M.  Acc. Chem. Res. 1972, 5, 225.
QUIVER: Saunders, M.; Laidig, K.E. Wolfsberg, M.  JACS 1989 111 8989.
Scaling Factors: Wong.  Chem. Phys. Lett. 1996, 256, 391-399.  Radom.  J Phys. Chem. 1996, 100, 16502.
Tunnelling Corrections: Bell.  Chem. Soc. Rev.  1974, 3, 513.

# Acknowledgements

The QUIVER program is modified from the original source code by Keith Laidig.
The scripts are refactored from the Houk group's qcrunch script and awk scripts from
Dan Singleton.

If you need help, please email me (ekwan16@gmail.com).

Eugene Kwan, November 2015
