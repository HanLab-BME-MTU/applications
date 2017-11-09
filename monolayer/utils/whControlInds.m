function [indsControl] = whControlInds(strLabels)
strPSup = 'pSuper';
strNT = 'NT';
strNT2 = 'N/T';
strDMSO = 'DMSO';
strDMSO24 = 'DMSO24';
strH2O = 'H2O';

indsControl = ...
    strcmp(whToPrefix(strLabels),strPSup) | strcmp(whToPrefix(strLabels),strNT) | strcmp(whToPrefix(strLabels),strNT2) | ...
    strcmp(whToPrefix(strLabels),strDMSO) | strcmp(whToPrefix(strLabels),strDMSO24) | strcmp(whToPrefix(strLabels),strH2O); 
end