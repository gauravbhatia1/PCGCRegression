# PCGCRegression
This is a low memory implementation of HE regression / PCGC Regression (Golan et al. 2014 PNAS)

It has 3 major changes relative to default for the purposes of performance:

1. The regression is performed on the fly meaning that O(k) memory is what is required for k variance components, instead of O(N^2).
2. PCs are subtracted from the GRM on the fly, saving the need to store multiple copies of the same GRM (pre- and post- PC subtraction)
3. Jackknifes are performed on the fly. This Jackknife is done over individuals and correction for fixed effects is not done with each jacknife. This increases the memory usage to O(k*l) for k variance components and l jackknife partitions.

Together, these changes make this method tractable for large, (N=50k) data-sets.

Usage:

java -jar PCGCRegression.jar parfile

parfile has the following parameters in key=value format:

- grmlist: Required. exactly the same file type as the "mgrm" file from the GCTA software. The individuals must be indentically ordered in all of the *.grm.id files.
- subtractvecs: Optional. a list of files corresponding to the files in grmlist. Each line must give a prefex for with <prefix>.eigenvec and <prefix>.eigenval exist. This must have exactly the same number of lines as grmlist.
- phenos: Required. Exactly the same phenotype file as used by GCTA. This currently expects case control data.
- covars: Optional. Covariates to include in the regression. These will be used in a logistic model to compute the constants to include in the regression. Inclusion of eigenvectors in the subtractvecs files does NOT include them here. They should be included.
- prevalence: Required. The prevalence of the phenotype (e.g. 0.01, or 0.001 for a disease with 1% or 0.1% prevalence, respectively)
- modcount: Optional. Given n, print an update to the console every n'th pair of GRM values. A reasonable value is 10^6.
phenos";
