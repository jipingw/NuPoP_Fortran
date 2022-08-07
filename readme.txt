1. About NuPoP_Fortran

NuPoP_Fortran is the Fortran codes to predict the nucleosome positioning based on DNA sequences. It is the core code from the NuPoP R package available from http://bioconductor.org/packages/release/bioc/html/NuPoP.html. This package provides a simple command line interface for nucleosome positioning prediction. More detailed information can be found in the NuPoP R package documentation enclosed (NuPoP.html).


2. Authors and maintainer

Author: Ji-Ping Wang <jzwang@northwestern.edu>; Liqun Xi <liqunxi02@gmail.com>
Maintainer: Ji-Ping Wang<jzwang@northwestern.edu>

3. License

License: GPL-2

3. License

License: GPL-2

4. How to?

a. Ensure you have Fortran compiler installed, e.g. generic fortran compiler (gfortran) or Intel fortran 
b. Compile the npred.f90 (mnase model) and npred_chem (chemimcal map model) from the command line, 
   
   $gfortran npred.f90 -o npred
   $gfortran npred_chem.f90 -o npred_chem
   
c. Run compiled Fortran codes:

   $ npred
or
   $ npred_chem

You will need to input sequence file name (in fasta format), species and etc. Two models are built in, 1st order duration HMM and 4th order duration HMM.
Results will typically be similar. The input of maximum linker length could 100 bp, 150 bp etc, which defines the possible maximum linker length. Longer length requries longer time, while the results will be similar provides a not-too-short linker length is specified. 100 bp should be appropriate for most applications.

5. Output

NuPoP_Fortran outputs five columns: chromosome coordinate (position), probability that a given chromosome location is the start of a nucleosome(P-start), nucleosome  occupancy(occup) at the given location, the nucleosome affinity score (N/L) and  histone affinity score for every 147 bp. npred will output file with extension _P1.txt  or P4.txt standing for first order and 4th order HMM. The npred_chem will output same files but with _chem extension.

6. References

- Wang, J.-P., Fondufe-Mittendorf, Y., Xi, L., Tsai, G.-F., Segal, E., and Widom, J. (2008). Prefer- entially quantized linker DNA lengths in Saccharomyces cerevisiae. PLoS Computational Biology, 4(9):e1000175.
 
- Xi, L., Fondufe-Mittendorf, Y., Xia, L., Flatow, J., Widom, J., and Wang, J.-P. (2010). Predict- ing nucleosome positioning using a duration hidden markov model. BMC Bioinformatics, pages doi:10.1186/1471–2105–11–346.


