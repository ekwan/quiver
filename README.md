## Scripts for Running QUIVER to Calculate KIEs

# Introduction

This is a series of scripts that facilitate the calculation of kinetic isotope effects using
Bigeleisen-Mayer theory.  You provide a configuration file containing the temperature of interest,
frequency scaling factor, and which isotopic replacements you want to make.  You then run the
enclosed scripts on two Gaussian output files containing the frequencies of the ground and
transition states  This script gives the KIEs, with and without standard tunnelling corrections.

# Installation Instructions
1. Clone this repository.  See [these instructions](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository).
2. Compile the QUIVER executable.
   * Run the `src/compile.sh` script.
   * This requires a Fortran 77 compiler.  The script assumes you are using gfortran. The provided
   executable was compiled on a MacBook.
   * Make sure the executable is called `quiver.exe` and is located in the same directory as the
   awk scripts.  This usually means the root directory of the repository.
3. Copy the demo files (`claisen_demo/*.out` and `claisen_demo/claisen.config`) to the root of the
repository.  Run the demonstration calculation with `./quiver.sh claisen.config claisen_gs.out claisen_ts.out`.
This assumes that you can run `bash` scripts.  This will work fine on Mac/Linux systems.  On Windows, use
Cygwin.

# Claisen Example

As a test case, I have provided the transition states for aliphatic Claisen rearrangement of
allyl vinyl ether that was studied by Singleton in *JACS* **1999**, *121*, 10865.  These files are
stored in the `claisen_demo/` folder.  There are three files.  The `claisen.config` file is the
input file for these scripts where the isotopic replacements are specified.  Two g09 `.out` files
are also provided, one for the ground state (gs) and one for the transition state (ts).  (Both are
calculated at B3LYP/6-31G(d).  I lifted the structures from the SI.  Conveniently, the carbons and
oxygen atom are exactly as labeled in Table 4 of the paper.)  The predicted values in the literature
for single atom <sup>12</sup>C/<sup>13</sup>C and <sup>16</sup>O/<sup>17</sup>O substitutions at 120
&deg;C are:

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
=== QUIVER Analysis ===

Temperature: 393 K
KIEs are referenced to isotopologue number 5.

Note: isotopomer descriptions refer to ground state atom numbers.

This is a kinetic isotope effect calculation.  Tunnelling corrections have been applied,
but may not be accurate for H/D/T KIEs.

isotopologue description                                  uncorrected      Widmer     infinite parabola
                                                              KIE           KIE              KIE
reference KIEs                                               1.002         1.002            1.002
isotopomer 1 - 13C @ atom 1                                  1.011         1.012            1.013
isotopomer 2 - 13C @ atom 2                                  1.000         1.000            1.000
isotopomer 3 - 17O @ atom 3                                  1.017         1.018            1.019
isotopomer 4 - 13C @ atom 4                                  1.028         1.031            1.031
[isotopomer 5 - 13C @ atom 5] (referenced to 1.00)           1.000         1.000            1.000
isotopomer 6 - 13C @ atom 6                                  1.013         1.015            1.015
isotopomer 7 - 2H @ atom 7, 2H @ atom 8                      0.953         0.954            0.955
```

As you can see, the numbers match to within 0.001.  (I'm not exactly sure which tunnelling
correction was used in Table 4.  Singleton notes that tunnelling improves KIE predictions
a bit.  For H/D predictions, though, these crude tunnelling corrections will be poor.  See
the Singleton Claisen paper for a fuller discussion, which references work by Truhlar in
this area.)  Note that this uses a scaling factor of 0.9614 for the frequencies.  Consult
the references section below for more scaling factors.

The exact input decks, temporary files, and expected output are also prvoided in the
`claisen_demo` folder.

# Gaussian Instructions

QUIVER needs frequencies from Gaussian to work.  This version has been tested with
Gaussian 09, revision D.  Optimize your geometries at your desired level of theory.
In a separate file (not strictly necessary), request a `#p freq` calculation at the same
level of theory.

* In theory, one can calculate KIEs at non-equilibrium geometries.  See
Singleton *JACS*  **2009**, *131*, 2397.  Note that the method used by QUIVER to calculate
frequencies is slightly different from that of g09, in that it doesn't project out the
translational and rotational modes.  Instead, it just ignores modes between -50 and 50 cm<sup>-1</sup>.
If running a grid, this could cause discontinuities as this script switches between EIE and KIE
calculations.  (The KIE calculations require the ratio of the imaginary frequencies if one structure
is a transition state.  This is like a classical correction for how quickly the system passes over
the barrier.)  The most rigorous thing to do would probably be to get the frequencies from g09
by performing frequency calculations from checkpoints and altering the isotopes with readisotopes.
Then, one could use the Bigeleisen-Mayer equation to calculate the KIE.  However, the work by
Singleton suggests this is unnecessary and this would be a very slow approach.
* Tighter convergence criteria can sometimes improve accuracy by reducing the size of
the smallest six translational/rotational frequencies.
* To calculate an equilibrium isotope effect, replace the transition state with the second ground state.
* The temperature you choose does not matter and will not be parsed from the file.

# Detailed Notes

Running the `quiver.sh` script starts the calculation:

`./quiver.sh claisen.config claisen_gs.out claisen_ts.out`

1. The script will first run the AWK script `quiver_prep.awk` to make the QUIVER input files.
The script first reads `claisen.config` and will then convert your instructions into the
format needed by QUIVER.  For instructions, see the comments in the file.  In theory, the
script will abort if you make a mistake, but this has not been tested perfectly.

2. QUIVER is run twice, once for the ground state and once for the transition state.
In each case, the mass-weighted Hessian is diagonalized repeatedly, first for the
unsubstituted isotopomer, and once for each desired substituted isotopomer.  The
Redlich-Teller product rule is then used to generate reduced isotopic partition
functions for each state.

3. The `quiver_analysis.awk` script analyzes the results.  This has only been tested
for KIEs.  The formula for the raw KIE is (imaginary_frequency<sub>light</sub> /
imaginary_frequency<sub>heavy</sub> x ( *f*(gs,heavy/gs,light) / *f*(ts,heavy/gs,light) ).
The imaginary frequency is just defined as the lowest (actual non-translational/non-vibrational)
vibrational frequency.  The analysis program just uses the imaginary frequency with the
largest magnitude for this ratio.  QUIVER ignores any negative frequencies, so large structures that
have small unconverged imgainary frequencies from numerical imprecision should work.  If there
are no negative frequencies, the analysis script will automatically perform an EIE calculation.
Note that the tunnelling corrections are not applicable in that case.

  *f* is a reduced isotopic partition function, which is marked as "(S2/S1)F" in
the QUIVER output.  gs and ts refer to the ground and transition states, respectively.

4. All temporary files are removed.  To prevent this behavior, comment out the
`rm` commands in the `quiver.sh` file.

# Technical Details

The scripts extract the last geometry, frequencies, Cartesian force constants, and some
other information from the two g09 .out files for the ground and transition state.  The
main Quiver input is dumped into `gs.q1` and `ts.q1`.  The Cartesian force constants are
placed in `gs.q2` and `ts.q2`.  Descriptions of the isotopologues being calculated are
placed in `gs.q3` and `ts.q3`.  These files are not used by QUIVER, but are needed for
the `quiver_analysis.awk` script.  Note that at present, only `gs.q3` is used (`ts.q3` is
ignored) to calculate the KIEs.

As much as possible, I designed the scripts to abort if any part of the parsing didn't
work or got data that make no sense.

# Compatability

This works on my Mac and has only been tested there.  (It is not computataionally demanding
to run--the expensive part of this is the frequencies, which you have to get from Gaussian).
It should work on any Linux system or Cygwin.  I have not tested it with the Windows version
of Gaussian.  The scripts depend on AWK, which should be available just about anywhere.
On Cygwin, be sure to install the GNU implementation of awk, gawk.  There could be an issue
running this with mawk because it does not seem to support left-justified printf.

# References

I lifted these references from the Houk group writeup for QUIVER.

1. **Bigeleisen-Mayer theory:**
  * J. Chem. Phys.*  **1947**, *15*, 261.
  * Wolfsberg, M.  *Acc. Chem. Res.* **1972**, *5*, 225.
2. **QUIVER:**
  * Saunders, M.; Laidig, K.E. Wolfsberg, M.  *JACS* **19898*, *111*, 8989.
3. **Scaling Factors:**
  * Wong.  *Chem. Phys. Lett.* **1996**, *256*, 391-399.
  * Radom.  *J Phys. Chem.* **1996**, *100*, 16502.
4. **Tunnelling Corrections:**
  * Bell.  *Chem. Soc. Rev.*  **1974**, *3*, 513.

# Troubleshooting

1. This program has been fixed to deal with up to 500 atoms and 100 isotopomers.  If there are further
issues, check the MNA and MMOL parameters in the quiver.f, as well as the input format string on line 92.

2. Please note that the atom numbers in the starting material and product do not need to be aligned, since
this calculates the reduced isotopic partition functions for the starting material and product separately.
The reasoning is as follows.  Suppose we have two reactions:

  A<sub>1</sub> + B --> A<sub>1</sub>B<sup>TS</sup>

  A<sub>2</sub> + B --> A<sub>2</sub>B<sup>TS</sup>

  where 1 = light and 2 = heavy.

  From transition state theory, the KIE for this reaction is K<sub>1</sub><sup>TS</sup>/K<sub>2</sub><sup>TS</sup>.
  From statistical mechanics, the equilibrium constant is given by the ratio of the partition functions:

  KIE = Q(A<sub>2</sub>)/Q(A<sub>1</sub>) x Q(B)/Q(B) x Q(A<sub>2</sub>B<sup>TS</sup>)/Q(A<sub>1</sub>B<sup>TS</sup>)

  Of course, the Q(B)/Q(B) term drops out.  QUIVER is run two batches, once to calculate the first ratio for each
  isotopomer, and again to calculate the third ratio for each isotopomer.  In brief, each partition function is
  assumed to be separable into translational, rotational, and vibrational components.  The Redlich-Teller rule is
  then used to construct a reduced partition function that is only terms of the component frequencies.  Each function
  is a product over each non-negative frequencies.  In the case of calculating a KIE, the ratio of the imaginary
  frequencies in the transition vectory is also required (multiplicatively).

  However, each isotope replacement must correspond to the same atom in the starting material and product.
  These atoms don't have to have the same number, but if you replace a hydrogen in one molecule and a carbon
  in the other, then the results will be gibberish.  This is checked in the `quiver_prep.awk` script.

3. If you are using any unusual atoms or isotopes, check the `quiver.f` and `quiver_prep.awk` files to
see if the necessary weights and symbols are present.

4. If there are strange errors when compiling QUIVER, check that the whitespace is exactly correct.

# Acknowledgements

The QUIVER program is modified from the original source code by Keith Laidig.
The scripts are refactored from the Houk group's qcrunch script and awk scripts from
Professor Daniel Singleton (Texas A&M).  I thank "Doc," as he is affectionately known,
for years of mentorship.

If you need help, please email me (ekwan16@gmail.com).

Eugene Kwan, December 2015
