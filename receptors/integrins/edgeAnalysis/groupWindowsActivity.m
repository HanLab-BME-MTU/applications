function sliceActivityGroup = groupWindowsActivity(protSamples,doPlot,...
    indxSlices,indxFrames,edgePosStd)
%GROUPWINDOWSACTIVITY groups windows based on activity categories
%
%SYNPOSIS sliceActivityGroup = groupWindowsActivity(protSamples,doPlot,...
%    indxSlices,indxFrames)
%
%INPUT  protSample : The protrusion samples as output by the windowing
%                    software.
%       doPlot     : Flag with value 1 to plot and 0 not to plot
%                    classification. Optional. Default: 0.
%       indxSlices : Vector with indices of windows (i.e. window number 
%                    along cell edge) to include in analysis.
%                    Optional. Default: all windows ([]).
%       indxFrames : Vector with indices of frames to include in analysis.
%                    Optional. Default: all frames ([]).
%       edgePosStd : Standard deviation of edge position.
%                    Optional. Default: 1.
%OUTPUT sliceActivityGroup: 27-by-1 structure array with field:
%            .edgeClassInfo: (number of events)-by-4 array. The 4
%                            columns store the frame index, slice index,
%                            length of event, and length of previous
%                            event.
%                    The 27 entries correspond to:
%                    Protrusion ...
%                    (1) Protrusion after retraction.
%                    (2) Protrusion after unknown.
%                    (3) Protrusion after short pause (1 frame), where pause is after protrusion.
%                    (4) Protrusion after short pause (1 frame), where pause is after retraction.
%                    (5) Protrusion after short pause (1 frame), where pause is after unknown.
%                    (6) Protrusion after long pause (>=2 frames), where pause is after protrusion.
%                    (7) Protrusion after long pause (>=2 frames), where pause is after retraction.
%                    (8) Protrusion after long pause (>=2 frames), where pause is after unknown.
%                    (9) Protrusion after anything, i.e. (1)-(8).
%                    Retraction ...
%                   (10) Retraction after protrusion.
%                   (11) Retraction after unknown.
%                   (12) Retraction after short pause (1 frame), where pause is after protrusion.
%                   (13) Retraction after short pause (1 frame), where pause is after retraction.
%                   (14) Retraction after short pause (1 frame), where pause is after unknown.
%                   (15) Retraction after long pause (>=2 frames), where pause is after protrusion.
%                   (16) Retraction after long pause (>=2 frames), where pause is after retraction.
%                   (17) Retraction after long pause (>=2 frames), where pause is after unknown.
%                   (18) Retraction after antyhing, i.e. (10)-(17).
%                    Pause ...
%                   (19) Pause after protrusion and before protrusion.
%                   (20) Pause after protrusion and before retraction.
%                   (21) Pause after protrusion and before unknown.
%                   (22) Pause after retraction and before protrusion.
%                   (23) Pause after retraction and before retraction.
%                   (24) Pause after retraction and before unknown.
%                   (25) Pause after unknown and before protrusion.
%                   (26) Pause after unknown and before retraction.
%                   (27) Pause after unknown and before unknown.
%
%Khuloud Jaqaman, January 2012

%% Input

if nargin < 1
    error('--groupWindowsActivity: Please enter protrusion samples!');
end

if nargin < 2 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 3 || isempty(indxSlices)
    indxSlices = [];
end

if nargin < 4 || isempty(indxFrames)
    indxFrames = [];
end

if nargin < 5 || isempty(edgePosStd)
    edgePosStd = 1;
end

%% Grouping

%classify the activity of each window at each interval
sliceMotionTypePause = classifyEdgeMotion(protSamples,doPlot,indxSlices,indxFrames,1,[],edgePosStd);
sliceMotionTypeProt = classifyEdgeMotion(protSamples,doPlot,indxSlices,indxFrames,2,[],edgePosStd);
sliceMotionTypeRetr = classifyEdgeMotion(protSamples,doPlot,indxSlices,indxFrames,3,[],edgePosStd);

%get number of window slices
numSlices = size(sliceMotionTypePause,1);

%initialize output structure
numTypes = 27;
sliceActivityGroup = repmat(struct('edgeClassInfo',NaN(0,4)),numTypes,1);

for iWindow = 1 : numSlices
    
    %extract this window's classification time series
    classCurrentPause = squeeze(sliceMotionTypePause(iWindow,:,:));
    classCurrentProt = squeeze(sliceMotionTypeProt(iWindow,:,:));
    classCurrentRetr = squeeze(sliceMotionTypeRetr(iWindow,:,:));
    
    %% protrusion
    
    %find protrusion onsets
    protOnset = find( classCurrentProt(:,1)==1 & diff([NaN; classCurrentProt(:,1)])~=0 );
    
    %keep only protrusions that last for >= 2 frames
    protOnset = protOnset(classCurrentProt(protOnset,4)>=2);
    
    %divide protrusion onsets based on the preceding event:
    protAftRetr = protOnset(classCurrentProt(protOnset,2)==2);
    protAftUnkn = protOnset(isnan(classCurrentProt(protOnset,2)));
    protAftShortPause = protOnset(classCurrentProt(protOnset,2)==3&classCurrentProt(protOnset,5)==1);
    protAftLongPause = protOnset(classCurrentProt(protOnset,2)==3&classCurrentProt(protOnset,5)>=2);
    
    %sub-divide protrusion onsets after pause further based on what
    %happened before the pause
    protAftShortPauseAftProt = protAftShortPause(classCurrentProt(protAftShortPause-1,2)==1);
    protAftShortPauseAftRetr = protAftShortPause(classCurrentProt(protAftShortPause-1,2)==2);
    protAftShortPauseAftUnkn = protAftShortPause(isnan(classCurrentProt(protAftShortPause-1,2)));
    protAftLongPauseAftProt  = protAftLongPause( classCurrentProt(protAftLongPause-1,2)==1);
    protAftLongPauseAftRetr  = protAftLongPause( classCurrentProt(protAftLongPause-1,2)==2);
    protAftLongPauseAftUnkn  = protAftLongPause( isnan(classCurrentProt(protAftLongPause-1,2)));
    
    %add to output structure
    if ~isempty(protAftRetr)
        sliceActivityGroup(1).edgeClassInfo = ...
            [sliceActivityGroup(1).edgeClassInfo; ...
            [protAftRetr iWindow*ones(length(protAftRetr),1) ...
            classCurrentProt(protAftRetr,4:5)]];
    end
    if ~isempty(protAftUnkn)
        sliceActivityGroup(2).edgeClassInfo = ...
            [sliceActivityGroup(2).edgeClassInfo; ...
            [protAftUnkn iWindow*ones(length(protAftUnkn),1) ...
            classCurrentProt(protAftUnkn,4) NaN(length(protAftUnkn),1)]];
    end
    if ~isempty(protAftShortPauseAftProt)
        sliceActivityGroup(3).edgeClassInfo = ...
            [sliceActivityGroup(3).edgeClassInfo; ...
            [protAftShortPauseAftProt iWindow*ones(length(protAftShortPauseAftProt),1)] ...
            classCurrentProt(protAftShortPauseAftProt,4:5)];
    end
    if ~isempty(protAftShortPauseAftRetr)
        sliceActivityGroup(4).edgeClassInfo = ...
            [sliceActivityGroup(4).edgeClassInfo; ...
            [protAftShortPauseAftRetr iWindow*ones(length(protAftShortPauseAftRetr),1)] ...
            classCurrentProt(protAftShortPauseAftRetr,4:5)];
    end
    if ~isempty(protAftShortPauseAftUnkn)
        sliceActivityGroup(5).edgeClassInfo = ...
            [sliceActivityGroup(5).edgeClassInfo; ...
            [protAftShortPauseAftUnkn iWindow*ones(length(protAftShortPauseAftUnkn),1)] ...
            classCurrentProt(protAftShortPauseAftUnkn,4:5)];
    end
    if ~isempty(protAftLongPauseAftProt)
        sliceActivityGroup(6).edgeClassInfo = ...
            [sliceActivityGroup(6).edgeClassInfo; ...
            [protAftLongPauseAftProt iWindow*ones(length(protAftLongPauseAftProt),1)]...
            classCurrentProt(protAftLongPauseAftProt,4:5)];
    end
    if ~isempty(protAftLongPauseAftRetr)
        sliceActivityGroup(7).edgeClassInfo = ...
            [sliceActivityGroup(7).edgeClassInfo; ...
            [protAftLongPauseAftRetr iWindow*ones(length(protAftLongPauseAftRetr),1)] ...
            classCurrentProt(protAftLongPauseAftRetr,4:5)];
    end
    if ~isempty(protAftLongPauseAftUnkn)
        sliceActivityGroup(8).edgeClassInfo = ...
            [sliceActivityGroup(8).edgeClassInfo; ...
            [protAftLongPauseAftUnkn iWindow*ones(length(protAftLongPauseAftUnkn),1)]...
            classCurrentProt(protAftLongPauseAftUnkn,4:5)];
    end
    
    %% retraction
    
    %find retraction onsets
    retrOnset = find( classCurrentRetr(:,1)==2 & diff([NaN; classCurrentRetr(:,1)])~=0 );
    
    %keep only retractions that last for >= 2 frames
    retrOnset = retrOnset(classCurrentRetr(retrOnset,4)>=2);
    
    %divide retraction onsets based on the preceding event:
    retrAftProt = retrOnset(classCurrentRetr(retrOnset,2)==1);
    retrAftUnkn = retrOnset(isnan(classCurrentRetr(retrOnset,2)));
    retrAftShortPause = retrOnset(classCurrentRetr(retrOnset,2)==3&classCurrentRetr(retrOnset,5)==1);
    retrAftLongPause = retrOnset(classCurrentRetr(retrOnset,2)==3&classCurrentRetr(retrOnset,5)>=2);
    
    %sub-divide retraction onsets after pause further based on what
    %happened before the pause
    retrAftShortPauseAftProt = retrAftShortPause(classCurrentRetr(retrAftShortPause-1,2)==1);
    retrAftShortPauseAftRetr = retrAftShortPause(classCurrentRetr(retrAftShortPause-1,2)==2);
    retrAftShortPauseAftUnkn = retrAftShortPause(isnan(classCurrentRetr(retrAftShortPause-1,2)));
    retrAftLongPauseAftProt  = retrAftLongPause( classCurrentRetr(retrAftLongPause-1,2)==1);
    retrAftLongPauseAftRetr  = retrAftLongPause( classCurrentRetr(retrAftLongPause-1,2)==2);
    retrAftLongPauseAftUnkn  = retrAftLongPause( isnan(classCurrentRetr(retrAftLongPause-1,2)));
    
    %add to output structure
    if ~isempty(retrAftProt)
        sliceActivityGroup(10).edgeClassInfo = ...
            [sliceActivityGroup(10).edgeClassInfo; ...
            [retrAftProt iWindow*ones(length(retrAftProt),1) ...
            classCurrentRetr(retrAftProt,4:5)]];
    end
    if ~isempty(retrAftUnkn)
        sliceActivityGroup(11).edgeClassInfo = ...
            [sliceActivityGroup(11).edgeClassInfo; ...
            [retrAftUnkn iWindow*ones(length(retrAftUnkn),1) ...
            classCurrentRetr(retrAftUnkn,4) NaN(length(retrAftUnkn),1)]];
    end
    if ~isempty(retrAftShortPauseAftProt)
        sliceActivityGroup(12).edgeClassInfo = ...
            [sliceActivityGroup(12).edgeClassInfo; ...
            [retrAftShortPauseAftProt iWindow*ones(length(retrAftShortPauseAftProt),1)] ...
            classCurrentRetr(retrAftShortPauseAftProt,4:5)];
    end
    if ~isempty(retrAftShortPauseAftRetr)
        sliceActivityGroup(13).edgeClassInfo = ...
            [sliceActivityGroup(13).edgeClassInfo; ...
            [retrAftShortPauseAftRetr iWindow*ones(length(retrAftShortPauseAftRetr),1)] ...
            classCurrentRetr(retrAftShortPauseAftRetr,4:5)];
    end
    if ~isempty(retrAftShortPauseAftUnkn)
        sliceActivityGroup(14).edgeClassInfo = ...
            [sliceActivityGroup(14).edgeClassInfo; ...
            [retrAftShortPauseAftUnkn iWindow*ones(length(retrAftShortPauseAftUnkn),1)] ...
            classCurrentRetr(retrAftShortPauseAftUnkn,4:5)];
    end
    if ~isempty(retrAftLongPauseAftProt)
        sliceActivityGroup(15).edgeClassInfo = ...
            [sliceActivityGroup(15).edgeClassInfo; ...
            [retrAftLongPauseAftProt iWindow*ones(length(retrAftLongPauseAftProt),1)]...
            classCurrentRetr(retrAftLongPauseAftProt,4:5)];
    end
    if ~isempty(retrAftLongPauseAftRetr)
        sliceActivityGroup(16).edgeClassInfo = ...
            [sliceActivityGroup(16).edgeClassInfo; ...
            [retrAftLongPauseAftRetr iWindow*ones(length(retrAftLongPauseAftRetr),1)] ...
            classCurrentRetr(retrAftLongPauseAftRetr,4:5)];
    end
    if ~isempty(retrAftLongPauseAftUnkn)
        sliceActivityGroup(17).edgeClassInfo = ...
            [sliceActivityGroup(17).edgeClassInfo; ...
            [retrAftLongPauseAftUnkn iWindow*ones(length(retrAftLongPauseAftUnkn),1)] ...
            classCurrentRetr(retrAftLongPauseAftUnkn,4:5)];
    end
    
    %% pause
    
    %find pause onsets
    pauseOnset = find( classCurrentPause(:,1)==3 & diff([NaN; classCurrentPause(:,1)])~=0 );
    
    %keep only pauses that last for >= 2 frames
    pauseOnset = pauseOnset(classCurrentPause(pauseOnset,4)>=2);
        
    %divide pause onsets based on events before and after
    pauseAftProtBefProt = pauseOnset(classCurrentPause(pauseOnset,2)==1&classCurrentPause(pauseOnset,3)==1);
    pauseAftProtBefRetr = pauseOnset(classCurrentPause(pauseOnset,2)==1&classCurrentPause(pauseOnset,3)==2);
    pauseAftProtBefUnkn = pauseOnset(classCurrentPause(pauseOnset,2)==1&isnan(classCurrentPause(pauseOnset,3)));
    pauseAftRetrBefProt = pauseOnset(classCurrentPause(pauseOnset,2)==2&classCurrentPause(pauseOnset,3)==1);
    pauseAftRetrBefRetr = pauseOnset(classCurrentPause(pauseOnset,2)==2&classCurrentPause(pauseOnset,3)==2);
    pauseAftRetrBefUnkn = pauseOnset(classCurrentPause(pauseOnset,2)==2&isnan(classCurrentPause(pauseOnset,3)));
    pauseAftUnknBefProt = pauseOnset(isnan(classCurrentPause(pauseOnset,2))&classCurrentPause(pauseOnset,3)==1);
    pauseAftUnknBefRetr = pauseOnset(isnan(classCurrentPause(pauseOnset,2))&classCurrentPause(pauseOnset,3)==2);
    pauseAftUnknBefUnkn = pauseOnset(isnan(classCurrentPause(pauseOnset,2))&isnan(classCurrentPause(pauseOnset,3)));
    
    %add to output structure
    if ~isempty(pauseAftProtBefProt)
        sliceActivityGroup(19).edgeClassInfo = ...
            [sliceActivityGroup(19).edgeClassInfo; ...
            [pauseAftProtBefProt iWindow*ones(length(pauseAftProtBefProt),1) ...
            classCurrentPause(pauseAftProtBefProt,4:5)]];
    end
    if ~isempty(pauseAftProtBefRetr)
        sliceActivityGroup(20).edgeClassInfo = ...
            [sliceActivityGroup(20).edgeClassInfo; ...
            [pauseAftProtBefRetr iWindow*ones(length(pauseAftProtBefRetr),1) ...
            classCurrentPause(pauseAftProtBefRetr,4:5)]];
    end
    if ~isempty(pauseAftProtBefUnkn)
        sliceActivityGroup(21).edgeClassInfo = ...
            [sliceActivityGroup(21).edgeClassInfo; ...
            [pauseAftProtBefUnkn iWindow*ones(length(pauseAftProtBefUnkn),1) ...
            classCurrentPause(pauseAftProtBefUnkn,4:5)]];
    end
    if ~isempty(pauseAftRetrBefProt)
        sliceActivityGroup(22).edgeClassInfo = ...
            [sliceActivityGroup(22).edgeClassInfo; ...
            [pauseAftRetrBefProt iWindow*ones(length(pauseAftRetrBefProt),1) ...
            classCurrentPause(pauseAftRetrBefProt,4:5)]];
    end
    if ~isempty(pauseAftRetrBefRetr)
        sliceActivityGroup(23).edgeClassInfo = ...
            [sliceActivityGroup(23).edgeClassInfo; ...
            [pauseAftRetrBefRetr iWindow*ones(length(pauseAftRetrBefRetr),1) ...
            classCurrentPause(pauseAftRetrBefRetr,4:5)]];
    end
    if ~isempty(pauseAftRetrBefUnkn)
        sliceActivityGroup(24).edgeClassInfo = ...
            [sliceActivityGroup(24).edgeClassInfo; ...
            [pauseAftRetrBefUnkn iWindow*ones(length(pauseAftRetrBefUnkn),1) ...
            classCurrentPause(pauseAftRetrBefUnkn,4:5)]];
    end
    if ~isempty(pauseAftUnknBefProt)
        sliceActivityGroup(25).edgeClassInfo = ...
            [sliceActivityGroup(25).edgeClassInfo; ...
            [pauseAftUnknBefProt iWindow*ones(length(pauseAftUnknBefProt),1) ...
            classCurrentPause(pauseAftUnknBefProt,4) NaN(length(pauseAftUnknBefProt),1)]];
    end
    if ~isempty(pauseAftUnknBefRetr)
        sliceActivityGroup(26).edgeClassInfo = ...
            [sliceActivityGroup(26).edgeClassInfo; ...
            [pauseAftUnknBefRetr iWindow*ones(length(pauseAftUnknBefRetr),1) ...
            classCurrentPause(pauseAftUnknBefRetr,4) NaN(length(pauseAftUnknBefRetr),1)]];
    end
    if ~isempty(pauseAftUnknBefUnkn)
        sliceActivityGroup(27).edgeClassInfo = ...
            [sliceActivityGroup(27).edgeClassInfo; ...
            [pauseAftUnknBefUnkn iWindow*ones(length(pauseAftUnknBefUnkn),1) ...
            classCurrentPause(pauseAftUnknBefUnkn,4) NaN(length(pauseAftUnknBefUnkn),1)]];
    end
    
end

%% combination groups
sliceActivityGroup(9).edgeClassInfo = vertcat(sliceActivityGroup(1:8).edgeClassInfo);
sliceActivityGroup(18).edgeClassInfo = vertcat(sliceActivityGroup(10:17).edgeClassInfo);

%% ~~~ the end ~~~

