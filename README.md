# PCGCRegression
This is a low memory implementation of HE regression / PCGC Regression (Golan et al. 2014 PNAS)

It has 3 major changes relative to default for the purposes of performance:

1. The regression is performed on the fly meaning that O(k) memory is what is required for k variance components, instead of O(N^2).
2. PCs are subtracted from the GRM on the fly, saving the need to store multiple copies of the same GRM (pre- and post- PC subtraction)
3. Jackknifes are performed on the fly. This Jackknife is done over individuals and correction for fixed effects is not done with each jacknife. This increases the memory usage to O(k*l) for k variance components and l jackknife partitions.

Together, these changes make this method tractable for large, (N=50k) data-sets.
