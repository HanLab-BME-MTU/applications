This contains Matlab functions for calculation of
distance transform of binary 3D images, bwdistsc 
and bwdistsc1. bwdistsc1 accepts additional parameter
to restrict the calculation of the distance transform 
only to a smaller band around the object, out to some 
distance maxval. This may significantly accelerate 
calculations in strongly non-convex situations, where 
the distance transform only out to a certain smaller 
value is required. 

Also included in this submission is the experimental 
version of bwdistsc, bwdistX, which uses an optimized 
forward-backward scan formulation of the algorithm, 
and is up to 5 times faster than bwdistsc (up to 10 
times faster than native bwdist). The function is 
currently included for user testing, and in the near 
future will replace the older bwdistsc in the package.
 
If you find any bugs or errors in the package, please 
report them to the author at the email address listed 
with Matlab Exchange.

Yuriy Mishchenko, PhD