function [filoInfo, filoWindowIdx ] = GCAassociateVeilWindowsToFilo(filoInfo,windows,plot)
% For each frame mark the windows with a filopodia so can exract protrusion
%/retract values specifically around the filo, record the window number in
%the filoInfo structure format. (was markWindowsWithFilo until 20140908)
% INPUT: 
% filoInfo: an N x 1 structure where N is the number of detected
% filopodia in a given frame and the structure fields represent a variety of
% information associated with each filopoda 
% windows: window as output by Hunter's windowing software for a given frame 
% 
% OUTPUT: 
% filoInfo:  adds a new field to the filoInfo structure .windowIdx 
%            a nx1 vector containing the slice idx of all corresponding
%            lamellipodia windows for that filo - where n is the number of 
%            windows within 1 pixel of the base of the filo coordinate 
% 
% filoWindowIdx: logical indexing marking which windows contain filopodia 

% xBase = arrayfun(@(x) x.Ext_coordsXY(1,1),filoInfo); 
% yBase = arrayfun(@(x) x.Ext_coordsXY(1,2),filoInfo);
% coordsBase = [xBase;yBase]; 
filoWindowIdx = zeros(numel(windows),1);
% maybe not most efficient way but works 
% for each filo find the nearest window. 
for ifilo = 1:numel(filoInfo) 
    coordsBase = [filoInfo(ifilo).Ext_coordsXY(1,1),filoInfo(ifilo).Ext_coordsXY(1,2)]; 
    coordsBase = coordsBase'; 
   % get a mat with windows of filopodia marked 
   iAssocWind = 1; % start count for windows- reset for each filo ( most will ony have one some might have multiple-decide later how want to deal) 
for iSlice = 1:numel(windows)-1
    % find the idx of windows with 1 pixel of the filo coord 
    if (~isempty(windows{iSlice}) &&  ~isempty(windows{iSlice}{1}))
        x = correspondingIndices(coordsBase, windows{iSlice}{1}{end},1);
        if isnan(x) ~=1
            filoInfo(ifilo).windowIdx(iAssocWind) = iSlice;
            filoWindowIdx(iSlice) = 1; % mark that the window contains a filo
            iAssocWind = iAssocWind+1;
            
        end
        clear x
    end
end
 
end 
% if plots plot filo specific windows 

if plot == 1
    filoWindowIdx = logical(filoWindowIdx); 
    plotWindows(windows(filoWindowIdx),{'r','FaceAlpha',.5},'bandMax',1); 
    plotWindows(windows(~filoWindowIdx),{'b','FaceAlpha',.5},'bandMax',1); 
    
end 


end

