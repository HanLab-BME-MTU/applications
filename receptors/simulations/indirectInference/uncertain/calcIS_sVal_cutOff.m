function [s,isVals_mod] = calcIS_sVal_cutOff(isVals_target,isVals_probe)
%CALCIS_SVAL_CUTOFF calculates the Mahalanobis distance S using the target
%and probe intermediate statistics given. This function is part of the 
%indirect inference based model calibration framework.
%
%   The current scheme is as follows-
%   Missing density values will be filled with 0s, but missing rates will
%   result in corresponding values in other vector/matrix to be removed.
%   Also not removing erroneous densities (densities for sizes > expected)
%   and using nancov instead of cov to accomodate the changes on
%   calcIS_cutOff.
%
%   INPUT:  
%           isVals_target:      a struct with target intermediate
%                               statistics and related values as output by
%                               calcIS_cutOff
%           isVals_probe:       a struct with probe intermediate
%                               statistics and related values as output by
%                               calcIS_cutOff
%
%   OUTPUT:
%           s:                  the calculated Mahalanobis distance
%           isVals_mod:         a struct with the following fields:
%                               1) theta_p - probe intermediate statistics
%                                            which may have been modified
%                               2) theta_t - target intermediate statistics
%                                            which may have been modified
%                               3) intStatsMatrix_p - probe intermediate 
%                                                     statistics matrix
%                                                     which may have been 
%                                                     modified
%                               4) intStatsMatrix_t - target intermediate 
%                                                     statistics matrix
%                                                     which may have been 
%                                                     modified
%                               5) intStatsMatrix_p - probe intermediate 
%                                                     statistics matrix
%                                                     which may have been 
%                                                     modified
%                               6) nancov_intStatsMatrix_t - covriance of 
%                                                     intStatsMatrix_t
%                               7) nancov_intStatsMatrix_p - covriance of 
%                                                     intStatsMatrix_p
%
%
%
%   Robel Yirdaw, 10/14/14.
%

    %Enable the following if using maximum cluster size restriction
    %The absolute largest cluster size expected, as set in the simulations
    %MAX_CLUST_SIZE = 5;
    
    theta_p = isVals_probe.theta;
    theta_t = isVals_target.theta;    
    
    intStatsMatrix_p = isVals_probe.intStatsMatrix;
    intStatsMatrix_t = isVals_target.intStatsMatrix;
    
    %The largestClustSize should not be greater than the absoule MAX
    %largestClustSize_p = min(isVals_probe.largestClustSize,MAX_CLUST_SIZE);
    %largestClustSize_t = min(isVals_target.largestClustSize,MAX_CLUST_SIZE);    

    %Get largest and cut-off cluster sizes for both target and probe
    largestClustSize_p = isVals_probe.largestClustSize;
    largestClustSize_t = isVals_target.largestClustSize;        
    
    cutOffClust_p = isVals_probe.cutOffClust;
    cutOffClust_t = isVals_target.cutOffClust;
    
    %The first set of intermediate statistics are densities. Missing values
    %will be filled with 0s.
    if (largestClustSize_t < largestClustSize_p)
        numMissingVals = largestClustSize_p - largestClustSize_t;
        %Reform theta
        theta_t = [theta_t(1,1:largestClustSize_t) zeros(1,numMissingVals)...
            theta_t(1,largestClustSize_t+1:end)];
        %Reform matrix
        intStatsMatrix_t = [intStatsMatrix_t(:,1:largestClustSize_t)...
            zeros(length(intStatsMatrix_t(:,1)),numMissingVals)...
            intStatsMatrix_t(:,largestClustSize_t+1:end)];        
        
    elseif (largestClustSize_p < largestClustSize_t)
        numMissingVals = largestClustSize_t - largestClustSize_p;
        %Reform theta
        theta_p = [theta_p(1,1:largestClustSize_p) zeros(1,numMissingVals)...
            theta_p(1,largestClustSize_p+1:end)];
        %Reform matrix
        intStatsMatrix_p = [intStatsMatrix_p(:,1:largestClustSize_p)...
            zeros(length(intStatsMatrix_p(:,1)),numMissingVals)...
            intStatsMatrix_p(:,largestClustSize_p+1:end)];                  
    end
    
    %Next set of intermediate statistics are the association and
    %dissociation rates
    if ( cutOffClust_t < cutOffClust_p )
              
        %Values in probe vector/matrix will be removed to match size of
        %target vector and matrix
        
        %Make a copy of the section to be modified (rates)
        tempISvec = theta_p(1,max(largestClustSize_t,largestClustSize_p)+1:end);
        %Reform theta
        theta_p = [theta_p(1,1:max(largestClustSize_t,largestClustSize_p))...
            tempISvec(1:cutOffClust_t) tempISvec(cutOffClust_p+1:cutOffClust_p+1+(cutOffClust_t - 2)) ];

        %Copy matrix
        tempISmat = intStatsMatrix_p(:,max(largestClustSize_t,largestClustSize_p)+1:end);        
        %Reform matrix
        intStatsMatrix_p = [intStatsMatrix_p(:,1:max(largestClustSize_t,largestClustSize_p) )...
            tempISmat(:,1:cutOffClust_t) tempISmat(:,cutOffClust_p+1:cutOffClust_p+1+(cutOffClust_t - 2)) ];
        
        
    elseif ( cutOffClust_p < cutOffClust_t )

        %Values in target vector/matrix will be removed to match size of
        %probe vector and matrix        

        %Make a copy of the section to be modified (rates)
        tempISvec = theta_t(1,max(largestClustSize_t,largestClustSize_p)+1:end);
        
        %Reform theta
        theta_t = [theta_t(1,1:max(largestClustSize_t,largestClustSize_p))...
            tempISvec(1:cutOffClust_p) tempISvec(cutOffClust_t+1:cutOffClust_t+1+(cutOffClust_p - 2)) ];

        %Copy matrix
        tempISmat = intStatsMatrix_t(:,max(largestClustSize_t,largestClustSize_p)+1:end);        
        %Reform matrix
        intStatsMatrix_t = [intStatsMatrix_t(:,1:max(largestClustSize_t,largestClustSize_p) )...
            tempISmat(:,1:cutOffClust_p) tempISmat(:,cutOffClust_t+1:cutOffClust_t+1+(cutOffClust_p - 2)) ];       
    end
            
    
    %Get covariance matrices
    vVec = {nancov(intStatsMatrix_p); nancov(intStatsMatrix_t)};
    
    %Calculate s
    s = (theta_p - theta_t)*(inv(vVec{1} + vVec{2}))*...
        transpose(theta_p - theta_t);
    
    %Return modified IS values.
    isVals_mod.theta_p = theta_p;
    isVals_mod.theta_t = theta_t;
    isVals_mod.intStatsMatrix_p = intStatsMatrix_p;
    isVals_mod.intStatsMatrix_t = intStatsMatrix_t;
    isVals_mod.nancov_intStatsMatrix_p = vVec{1};
    isVals_mod.nancov_intStatsMatrix_t = vVec{2};
        
    
end %function

        
        