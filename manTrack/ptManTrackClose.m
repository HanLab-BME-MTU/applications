function ptManTrackClose(src,eventdata)
% ptManTrackClose is the callback for the clickable cell nuclei
%
% SYNOPSIS       ptManTrackClose
%
% INPUT          none
%
% OUTPUT         none
%
% DEPENDENCIES   ptManTrackClose uses {nothing}
%                                  
%                ptManTrackClose is used by { ptShowImageSequence, ptShowNextInSequence }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Feb 05          First Release

% get a handle to the current axes
h = get(gca);

% Get a handle to the figure
h2 = get(gcf);

hFig = findall(0,'Tag','ptSequenceFigure');
if isempty(hFig)
    % Nothing to delete
    return;
else
    % Get a handle to the figure
    handles = guidata(hFig);
    
    % Get the handles struct of the main gui
    handlesMain = guidata(gcbo);

    % Get the current frame number
    imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'imageCounter');
    lastTrackedFrame = str2num(get(imageCounterHandle, 'String'));
    
    % Fetch the values from the handlesMt
    resultsDirectory = handlesMain.resultsDirectory;
    saveDirectory = handlesMain.saveDirectory;
    imageRange = handlesMain.imageRange;
    MPM = handlesMain.MPM;
    manualMPM = handles.manualMPM;
    validFrames = handlesMain.validFrames;
    bodyName = handlesMain.bodyName;
    
    % Check whether the user has already clicked on some cells. If not
    % there is really nothing to do and save
    if isfield(handlesMain,'firstTrackedFrame')
        
        firstTrackedFrame = handlesMain.firstTrackedFrame;

        % Initialize index array
        indxInMPM = zeros(2,length(validFrames));

        % Get the index of the cell we first started tracking
        [value,indx] = min(sqrt((MPM(:,2*firstTrackedFrame-1)-manualMPM(2*firstTrackedFrame-1)).^2 + ...
                                (MPM(:,2*firstTrackedFrame)-manualMPM(2*firstTrackedFrame)).^2));
        MPMindx = indx;
        
        %for frameCount = firstTrackedFrame : length(validFrames)
        for frameCount = 1 : length(validFrames)     

            % Find the coordinates of the first manually tracked coordinates in the
            % MPM matrix
            MPMCoords = MPM(:,2*frameCount-1:2*frameCount);

            % Get manual coords
            manCoords = manualMPM(2*frameCount-1:2*frameCount);

            if ~isempty(find(manCoords))
                % Pick out the ones that are closed to the clicked point
                [value,indx] = min(sqrt((MPMCoords(:,1)-manCoords(1)).^2 + ...
                                        (MPMCoords(:,2)-manCoords(2)).^2));

                % Store the index
                indxInMPM(1,frameCount) = indx;    
            else
                indxInMPM(1,frameCount) = 0;
            end
            
            if ~isempty(find(MPMCoords(MPMindx,1) & MPMCoords(MPMindx,2)))
                indxInMPM(2,frameCount) = MPMindx;
            else
                indxInMPM(2,frameCount) = 0;
            end
        end

        if MPMindx ~= 0
            disp(indxInMPM);
            save([saveDirectory filesep 'manTrackResults_' bodyName(1:end-2) '_cell' num2str(MPMindx) '.mat'], 'indxInMPM');
        end
    end
    
    % Destroy the General Properties GUI
    handles = guidata(hFig);
    delete(handles.ptSequenceFigure);
end



