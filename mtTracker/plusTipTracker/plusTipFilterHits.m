function [ hits ] = plusTipFilterHits( saveDir, nameOutput, hits, stringency)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<1 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for hit files.');
end

paramNames = fieldnames(hits); 

for iParam = 1:length(paramNames)
    
    
    hitParam =  hits.(char(paramNames(iParam)))(2:end,:) ; 
    pValues = cell2mat(hitParam(:,4)); 
    idx = pValues < stringency;
    
  hitLessStr = hitParam(idx,:); 
  
    
   if isempty(hitLessStr) ~= 1
hits.(char(paramNames(iParam))) = [hits.(char(paramNames(iParam)))(1,:); hitLessStr];
   else 
     hits = rmfield(hits,char(paramNames(iParam))); 
   end 
end 
save([saveDir filesep nameOutput],'hits'); 
end

