function [projList] = getLabelName(dataSetArray,projList,colorGroup)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



projListOrig = projList; 
clear projList; 
clickpoints = true; 
count = 1; 
idxSave = []; 
while clickpoints == true
reply2 = questdlg('Document Points ?');


if strcmpi(reply2,'yes')
    
    hpoint = impoint;
    p = wait(hpoint);
    %p = round(p);
    
    
    % get the axis labels
    hLabelX = get(gca,'xLabel');
    nameX = get(hLabelX,'String');
    hLabelY = get(gca,'yLabel');
    nameY  = get(hLabelY,'String');
    
    % get the xvalue
    varName = get(dataSetArray,'VarNames');
    nameX = strrep(nameX,' ','_');
    
    nameX =  strrep(nameX,'filo','F');
    nameX = strrep(nameX,'Length','L');
    nameX =  strrep(nameX,'Intensity','I');
    nameX = strrep(nameX,'retractionAnalysis_perTime','RPers');
    nameX = strrep(nameX,'protrusionAnalysis_perTime','PPers');
    nameX = strrep(nameX,'protrusionAnalysis_mednVel','PVel');
    nameX =  strrep(nameX,'retractionAnalysis_mednVel','RVel');
    
    
    
    
    c1 = find(strcmpi(nameX,varName));
    c2 = numel(varName); %for now just set it to end
    
    
    
    % find the ones closest to
    xValues = dataSetArray.(varName{c1});
    yValues = dataSetArray.(varName{c2});
    
    delt1 = xValues-p(1);
    delt2 = yValues-p(2) ;
    totDelt = abs(delt1)+abs(delt2);
    
    
    
    
    
    idxProj = totDelt == min(totDelt);
    nameProjs = get(dataSetArray,'ObsNames');
    nameProjSelect = nameProjs{idxProj};
    
    nameProjSelect = strrep(nameProjSelect,'_',' ');
    
    text(xValues(idxProj),yValues(idxProj),nameProjSelect);
   
    scatter(xValues(idxProj),yValues(idxProj),50,colorGroup,'filled'); 
    
    idxSave(count) = find(idxProj);
    
    count = count+1; 
    
elseif strcmpi(reply2,'no')  || strcmpi(reply2,'cancel')
    
    clickpoints = false; 
     
end
if ~isempty(idxSave)
projList = projListOrig(idxSave,:); 
end 
end % while 











end

