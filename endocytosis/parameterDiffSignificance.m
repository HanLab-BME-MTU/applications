function [kspval] = parameterDiffSignificance(merCTRL,merTR)
% determine significance in the difference of each fit parameter between
% treatment and control, using a KS-test
% SYNOPSIS: [kspval] = parameterDiffSignificance(merCTRL,merTR)
% INPUT     merCTRL =   merged control data
%           merTR =     merged treatment data
%           The merged data have to contain fields called 
%           .hvec_cum, .tvec_cum, and .compactFitRes
%           NOTE: The merged data structures are the result of putting the
%           original data structure through the analysis functions:
%           sData = stageFitLifetimesPlat(OriginalData);
%           merData = mergeFastSlowHistogramsPlat(sData, 300, [2 2 1]);
%
% OUTPUT    kspval =    p-value of the ks-test, for all parameters; this
%                       vector has the same length as the parameter vector
% 
% last modified: Dinah Loerke 05/20/2008


% extract fit parameters of the control results and write them into the
% vector parCTRL
pmatCTRL = merCTRL.compactFitRes;
parCTRL([1,2,5,8]) = 0.01*pmatCTRL(:,1);
parCTRL([3:3:9]) = pmatCTRL(2:4,2);
parCTRL([4:3:10]) = [2 2 1];

chistCTRL = merCTRL.hvec_cum;
tvecCTRL = merCTRL.tvec_cum;



% extract fit parameters of the treatment results and write them into the
% vector parTR
pmatTR = merTR.compactFitRes;
parTR([1,2,5,8]) = 0.01*pmatTR(:,1);
parTR([3:3:9]) = pmatTR(2:4,2);
parTR([4:3:10]) = [2 2 1];

chistTR = merTR.hvec_cum;
tvecTR = merTR.tvec_cum;

% prepare fixvector (the vector which determines which parameters are kept
% fixed during the fit)
gfixvec = 0*parTR;
gfixvec(4:3:10) = 1;


% loop over all fit parameters
for p=1:length(parTR)
    
    % fixed fit: set the currently selected parameter (e.g. tau1) to the
    % value of the control data, and fit the treatment data with this value
    % kept fixed
    startvec = parTR;
    startvec(p) = parCTRL(p);
    fixvec = gfixvec;
    fixvec(p) = 1;
    
    [estFix, resFix] = fitcurveMultiWeibullCDF_lsq(tvecTR, chistTR, startvec, fixvec);
    
    % free fit: using the results of the fixed fit as start parameters, now
    % fit again with allowing this parameter to converge freely
    [estFree, resFree] = fitcurveMultiWeibullCDF_lsq(tvecTR, chistTR, estFix,gfixvec);
    
    % compare the residuals of the free and the fixed fit to determine
    % whether the goodness-of-fit improved significantly by allowing this
    % parameter to deviate from the control value
    [h, kspval(p)] = kstest2(resFix,resFree);
    
    
   
    %%=====================================================================
    % special case (only for expert user):
    % if you are fitting only 2 populations in the treatment case, as 
    % opposed to 3 populations in the control (as is the case for cargo
    % OX), the uncomment the following paragraph, and read out the new
    % p-values by setting a break point at the end of the function or
    % overwriting kspval
    
    
%     startvec2 = startvec;
%     startvec2(5:7) = [];
%     
%     fixvec2 = fixvec;
%     fixvec2(5:7) = [];
%     gfixvec2 = gfixvec; 
%     gfixvec2(5:7) = [];
%     
%     [estFix2, resFix2] = fitcurveMultiWeibullCDF_lsq(tvecTR, chistTR, startvec2, fixvec2);
%     [estFree2, resFree2] = fitcurveMultiWeibullCDF_lsq(tvecTR, chistTR, estFix2,gfixvec2);
%     
%     [h, kspval2p(p)] = kstest2(resFix2,resFree2);
%     
    
    % end of special case
    %%=====================================================================
    
end


end % of function

    
    