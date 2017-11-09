function [ analInfo ] = addLocalVeilActivityField(rawData,analInfo) 
%small helper function for now to just add the local protrusion retraction to 
% the filoInfo 

for iFrame = 1:numel(analInfo)-1
    filoInfo = analInfo(iFrame).filoInfo; 
   


% add these local values to the analInfo 
for iFilo = 1:length(filoInfo)
    % get window 
    
    if ~isempty(filoInfo(iFilo).windowIdx)
        windNum = filoInfo(iFilo).windowIdx(1) ; 
    localVeil = rawData(windNum,iFrame); 
    
    filoInfo(iFilo).localVeil = localVeil; 
    else 
        filoInfo(iFilo).localVeil = NaN;
    end 
    
    
    
    
end 

analInfo(iFrame).filoInfo = filoInfo; 





end 
   
    
    
end 

