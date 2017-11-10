function isControl = isControlStr2016(strLabel)
strPSup = 'pSuper';
strNT = 'NT';
strDMSO = 'DMSO';
isControl = strcmp(strLabel,strPSup) | strcmp(strLabel,strNT) | strcmp(strLabel,strDMSO);
end