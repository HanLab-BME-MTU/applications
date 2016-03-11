function [ ratioValues,filoSetProt] = plotFiloActivityProtRetract(MD,analInfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 veilType{1} = 'retractionpersVeil'; 
 veilType{2} = 'protrusionpersVeil'; 
files = searchFiles('.mat',[],[MD.outputDirectory_ filesep 'ratio_images' filesep 'ratio_of_2_to_1'],0,'all',1); 
 
countProt = 1; 
for iFrame =1:numel(analInfo) 
    load(files{iFrame});
    % get the filo for current frame 
    filoInfo =  analInfo(iFrame).filoInfo; 
    
   for iVeilType = 1:numel(veilType)
    % get the indices of the filo 
    
    idxRetract = arrayfun(@(x) x.(veilType{iVeilType}),filoInfo); 
    idxLong = arrayfun(@(x) x.Ext_length>5,filoInfo); 
    idxSet = find(idxRetract& idxLong); 
    if ~isempty(idxSet)
    filoInfoRetract = filoInfo(idxSet);
    if (iVeilType == 1 && ~isempty(idxSet))
    for iFilo = 1:length(idxSet)
       filoSetProt(countProt).frameNum = iFrame; 
       filoSetProt(countProt).idxFilo = idxSet(iFilo); 
        countProt = countProt+1; 
       
    end
    end 
    end
    for i = 1:numel(filoInfoRetract)
        pixIndices = filoInfoRetract(i).('Ext_pixIndices');
        idxEnd = find(pixIndices == filoInfoRetract(i).('Ext_endpointCoordFitPix'));
        
        indicesFit{i} = pixIndices(1:idxEnd); 
       
    end
% load the FRET data 
    ratioValues{iVeilType}{iFrame} = cellfun(@(x) currRatio(x),indicesFit,'uniformoutput',0); 
    % if plot 
    
   end 
end 
plot =1 ; 

if plot == 1
    color(1) = 'r'; 
    color(2) = 'b'; 
    for iVeilType = 1:2
    values = horzcat(ratioValues{iVeilType}{:})';
    meanValues = cellfun(@(x) mean(x),values); 
    hist(meanValues,100); 
    hold on 
    end 
end 
% plot all mean value per filopodia from retracting and protruding regions.
% ratioValues{:}{1}

