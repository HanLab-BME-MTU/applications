function [xCenterOverTime,yCenterOverTime,areaOverTime] = ...
    windowCharacteristics(winPositions)
%WINDOWCHARACTERISTICS calculates window area and center over time
%
%SYNOPSIS [xCenterOverTime,yCenterOverTime,areaOverTime] = ...
%    windowCharacteristics(winPositions)
%
%INPUT  winPositions   : A 2D array of the window edges. 
%                        Number of rows = number of window frames. 
%                        Number of columns = number of windows parallel to
%                        the edge.
%                        Each entry is the output of Hunter's new
%                        windowing software.
%                        To make this variable, one puts together the
%                        windows of each frame coming out of the windowing
%                        software.
%
%OUTPUT xCenterOverTime: 3D array of dimensions (number of bands) x
%                        (number of slices) x (number of window frames)
%                        storing for each window in each frame the
%                        x-coordinate of its center.
%       yCenterOverTime: 3D array of dimensions (number of bands) x
%                        (number of slices) x (number of window frames)
%                        storing for each window in each frame the
%                        y-coordinate of its center.
%       areaOverTime   : 3D array of dimensions (number of bands) x
%                        (number of slices) x (number of window frames)
%                        storing for each window in each frame its area.
%
%Khuloud Jaqaman, August 2011

%% Input

if nargin < 1
    disp('--windowCharacteristics: Incorrect number of input arguments!');
    return
end

%% Pre-processing

%get number of frames that have windows and number of windows parallel to
%the edge
[numWinFrames,numWinPara] = size(winPositions);

%find number of windows perpendicular to the edge
nBands = cellfun(@(x)(numel(x)),winPositions);
numWinPerp = max(nBands(:));

%% Track assignment into windows

%initialize arrays of center coordinates and area
xCenterOverTime = NaN(numWinPerp,numWinPara,numWinFrames);
yCenterOverTime = NaN(numWinPerp,numWinPara,numWinFrames);
areaOverTime = NaN(numWinPerp,numWinPara,numWinFrames);

%go over all window frames
for iWinFrame = 1 : numWinFrames
    
    %go over the windows in this frame
    for iPara = 1 : numWinPara
        for iPerp = 1 : nBands(iWinFrame,iPara)
            
            %if this window has a finite size
            if ~isempty(winPositions{iWinFrame,iPara}{iPerp})
                
                %get the window boundaries
                windowsPoly = [winPositions{iWinFrame,iPara}{iPerp}{:}];
                winX = windowsPoly(1,:);
                winY = windowsPoly(2,:);
                
                %get window center coordinates
                xCenterOverTime(iPerp,iPara,iWinFrame) = mean(winX);
                yCenterOverTime(iPerp,iPara,iWinFrame) = mean(winY);
                
                %get window area
                areaOverTime(iPerp,iPara,iWinFrame) = polyarea(winX,winY);
                
            else %if this window is collapsed, then there are no tracks in it
                
                xCenterOverTime(iPerp,iPara,iWinFrame) = NaN;
                yCenterOverTime(iPerp,iPara,iWinFrame) = NaN;
                areaOverTime(iPerp,iPara,iWinFrame) = 0;
                
            end
            
        end
    end
    
end
