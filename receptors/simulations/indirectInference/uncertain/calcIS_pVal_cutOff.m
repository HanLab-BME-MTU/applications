function [p] = calcIS_pVal_cutOff(s,dof)
%CALCIS_PVAL_CUTOFF calculates p-value corresponding to the Mahalanobis
%distance s, using degree of freedom dof.  This function is called by
%getIS_spVals and is part of the indirect inference based model calibration
%framework.
%
%   INPUT:  
%           s:      Mahalanobis distance for target and probe IS
%           dof:    degree of freedom
%
%   OUTPUT:
%           p:      calculated p-value.
%
%   Robel Yirdaw, November 2014
%


    p = 1 - chi2cdf(s,dof);
    
end