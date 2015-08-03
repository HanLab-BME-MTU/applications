function [ scale ] = qqscale( A,B )
%Obtain scale factor similar to qqplot

scale = lsqr(prctile(A,25:75)',prctile(B,25:75)');


end

