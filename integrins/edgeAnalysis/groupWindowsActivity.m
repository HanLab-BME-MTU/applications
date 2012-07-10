function windowGroup = groupWindowsActivity(protSamples,doPlot,...
    windowRange,frameRange)
%GROUPWINDOWSACTIVITY groups windows based on activity categories
%
%SYNPOSIS windowGroup = groupWindowsActivity(protSamples,doPlot,...
%    windowRange,frameRange)
%
%INPUT  protSample : The protrusion samples as output by the windowing
%                    software.
%       doPlot     : Flag with value 1 to plot and 0 not to plot
%                    classification. Optional. Default: 0.
%       windowRange: 2-row vector indicating range of windows (i.e.
%                    window number along cell edge) to include in analysis.
%                    Optional. Default: [].
%       frameRange : 2-row vector indicating range of window frames
%                    to include in analysis.
%                    Optional. Default: [].
%OUTPUT windowGroup: 27-by-1 structure array with field:
%            .edgeClassInfo: (number of events)-by-4 array. The 4
%                            columns store the frame index, window index,
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

if nargin < 3 || isempty(windowRange)
    windowRange = [];
end

if nargin < 4 || isempty(frameRange)
    frameRange = [];
end

%% Grouping

%classify the activity of each window at each interval
windowMotionType = classifyEdgeMotion(protSamples,doPlot,windowRange,frameRange);

%get number of windows and frames
numWindows = size(windowMotionType,1);

%initialize output structure
numTypes = 27;
windowGroup = repmat(struct('edgeClassInfo',NaN(0,4)),numTypes,1);

for iWindow = 1 : numWindows
    
    %extract this window's classification time series
    classCurrent = squeeze(windowMotionType(iWindow,:,:));
    
    %% protrusion
    
    %find protrusion onsets
    protOnset = find( classCurrent(:,1)==1 & diff([NaN; classCurrent(:,1)])~=0 );
    
    %keep only protrusions that last for >= 2 frames
    protOnset = protOnset(classCurrent(protOnset,4)>=2);
    
    %divide protrusion onsets based on the preceding event:
    protAftRetr = protOnset(classCurrent(protOnset,2)==2);
    protAftUnkn = protOnset(isnan(classCurrent(protOnset,2)));
    protAftShortPause = protOnset(classCurrent(protOnset,2)==3&classCurrent(protOnset,5)==1);
    protAftLongPause = protOnset(classCurrent(protOnset,2)==3&classCurrent(protOnset,5)>=2);
    
    %sub-divide protrusion onsets after pause further based on what
    %happened before the pause
    protAftShortPauseAftProt = protAftShortPause(classCurrent(protAftShortPause-1,2)==1);
    protAftShortPauseAftRetr = protAftShortPause(classCurrent(protAftShortPause-1,2)==2);
    protAftShortPauseAftUnkn = protAftShortPause(isnan(classCurrent(protAftShortPause-1,2)));
    protAftLongPauseAftProt  = protAftLongPause( classCurrent(protAftLongPause-1,2)==1);
    protAftLongPauseAftRetr  = protAftLongPause( classCurrent(protAftLongPause-1,2)==2);
    protAftLongPauseAftUnkn  = protAftLongPause( isnan(classCurrent(protAftLongPause-1,2)));
    
    %add to output structure
    windowGroup(1).edgeClassInfo = ...
        [windowGroup(1).edgeClassInfo; ...
        [protAftRetr iWindow*ones(length(protAftRetr),1) ...
        classCurrent(protAftRetr,4:5)]];
    windowGroup(2).edgeClassInfo = ...
        [windowGroup(2).edgeClassInfo; ...
        [protAftUnkn iWindow*ones(length(protAftUnkn),1) ...
        classCurrent(protAftUnkn,4:5)]];
    windowGroup(3).edgeClassInfo = ...
        [windowGroup(3).edgeClassInfo; ...
        [protAftShortPauseAftProt iWindow*ones(length(protAftShortPauseAftProt),1)] ...
        classCurrent(protAftShortPauseAftProt,4:5)];
    windowGroup(4).edgeClassInfo = ...
        [windowGroup(4).edgeClassInfo; ...
        [protAftShortPauseAftRetr iWindow*ones(length(protAftShortPauseAftRetr),1)] ...
        classCurrent(protAftShortPauseAftRetr,4:5)];
    windowGroup(5).edgeClassInfo = ...
        [windowGroup(5).edgeClassInfo; ...
        [protAftShortPauseAftUnkn iWindow*ones(length(protAftShortPauseAftUnkn),1)] ...
        classCurrent(protAftShortPauseAftUnkn,4:5)];
    windowGroup(6).edgeClassInfo = ...
        [windowGroup(6).edgeClassInfo; ...
        [protAftLongPauseAftProt iWindow*ones(length(protAftLongPauseAftProt),1)]...
        classCurrent(protAftLongPauseAftProt,4:5)];
    windowGroup(7).edgeClassInfo = ...
        [windowGroup(7).edgeClassInfo; ...
        [protAftLongPauseAftRetr iWindow*ones(length(protAftLongPauseAftRetr),1)] ...
        classCurrent(protAftLongPauseAftRetr,4:5)];
    windowGroup(8).edgeClassInfo = ...
        [windowGroup(8).edgeClassInfo; ...
        [protAftLongPauseAftUnkn iWindow*ones(length(protAftLongPauseAftUnkn),1)]...
        classCurrent(protAftLongPauseAftUnkn,4:5)];
    
    %% retraction
    
    %find retraction onsets
    retrOnset = find( classCurrent(:,1)==2 & diff([NaN; classCurrent(:,1)])~=0 );
    
    %keep only retractions that last for >= 2 frames
    retrOnset = retrOnset(classCurrent(retrOnset,4)>=2);
    
    %divide retraction onsets based on the preceding event:
    retrAftProt = retrOnset(classCurrent(retrOnset,2)==1);
    retrAftUnkn = retrOnset(isnan(classCurrent(retrOnset,2)));
    retrAftShortPause = retrOnset(classCurrent(retrOnset,2)==3&classCurrent(retrOnset,5)==1);
    retrAftLongPause = retrOnset(classCurrent(retrOnset,2)==3&classCurrent(retrOnset,5)>=2);
    
    %sub-divide retraction onsets after pause further based on what
    %happened before the pause
    retrAftShortPauseAftProt = retrAftShortPause(classCurrent(retrAftShortPause-1,2)==1);
    retrAftShortPauseAftRetr = retrAftShortPause(classCurrent(retrAftShortPause-1,2)==2);
    retrAftShortPauseAftUnkn = retrAftShortPause(isnan(classCurrent(retrAftShortPause-1,2)));
    retrAftLongPauseAftProt  = retrAftLongPause( classCurrent(retrAftLongPause-1,2)==1);
    retrAftLongPauseAftRetr  = retrAftLongPause( classCurrent(retrAftLongPause-1,2)==2);
    retrAftLongPauseAftUnkn  = retrAftLongPause( isnan(classCurrent(retrAftLongPause-1,2)));
    
    %add to output structure
    windowGroup(10).edgeClassInfo = ...
        [windowGroup(10).edgeClassInfo; ...
        [retrAftProt iWindow*ones(length(retrAftProt),1) ...
        classCurrent(retrAftProt,4:5)]];
    windowGroup(11).edgeClassInfo = ...
        [windowGroup(11).edgeClassInfo; ...
        [retrAftUnkn iWindow*ones(length(retrAftUnkn),1) ...
        classCurrent(retrAftUnkn,4:5)]];
    windowGroup(12).edgeClassInfo = ...
        [windowGroup(12).edgeClassInfo; ...
        [retrAftShortPauseAftProt iWindow*ones(length(retrAftShortPauseAftProt),1)] ...
        classCurrent(retrAftShortPauseAftProt,4:5)];
    windowGroup(13).edgeClassInfo = ...
        [windowGroup(13).edgeClassInfo; ...
        [retrAftShortPauseAftRetr iWindow*ones(length(retrAftShortPauseAftRetr),1)] ...
        classCurrent(retrAftShortPauseAftRetr,4:5)];
    windowGroup(14).edgeClassInfo = ...
        [windowGroup(14).edgeClassInfo; ...
        [retrAftShortPauseAftUnkn iWindow*ones(length(retrAftShortPauseAftUnkn),1)] ...
        classCurrent(retrAftShortPauseAftUnkn,4:5)];
    windowGroup(15).edgeClassInfo = ...
        [windowGroup(15).edgeClassInfo; ...
        [retrAftLongPauseAftProt iWindow*ones(length(retrAftLongPauseAftProt),1)]...
        classCurrent(retrAftLongPauseAftProt,4:5)];
    windowGroup(16).edgeClassInfo = ...
        [windowGroup(16).edgeClassInfo; ...
        [retrAftLongPauseAftRetr iWindow*ones(length(retrAftLongPauseAftRetr),1)] ...
        classCurrent(retrAftLongPauseAftRetr,4:5)];
    windowGroup(17).edgeClassInfo = ...
        [windowGroup(17).edgeClassInfo; ...
        [retrAftLongPauseAftUnkn iWindow*ones(length(retrAftLongPauseAftUnkn),1)] ...
        classCurrent(retrAftLongPauseAftUnkn,4:5)];
    
    %% pause
    
    %find pause onsets
    pauseOnset = find( classCurrent(:,1)==3 & diff([NaN; classCurrent(:,1)])~=0 );
    
    %keep only pauses that last for >= 2 frames
    pauseOnset = pauseOnset(classCurrent(pauseOnset,4)>=2);
    
    %divide pause onsets based on events before and after
    pauseAftProtBefProt = pauseOnset(classCurrent(pauseOnset,2)==1&classCurrent(pauseOnset,3)==1);
    pauseAftProtBefRetr = pauseOnset(classCurrent(pauseOnset,2)==1&classCurrent(pauseOnset,3)==2);
    pauseAftProtBefUnkn = pauseOnset(classCurrent(pauseOnset,2)==1&isnan(classCurrent(pauseOnset,3)));
    pauseAftRetrBefProt = pauseOnset(classCurrent(pauseOnset,2)==2&classCurrent(pauseOnset,3)==1);
    pauseAftRetrBefRetr = pauseOnset(classCurrent(pauseOnset,2)==2&classCurrent(pauseOnset,3)==2);
    pauseAftRetrBefUnkn = pauseOnset(classCurrent(pauseOnset,2)==2&isnan(classCurrent(pauseOnset,3)));
    pauseAftUnknBefProt = pauseOnset(isnan(classCurrent(pauseOnset,2))&classCurrent(pauseOnset,3)==1);
    pauseAftUnknBefRetr = pauseOnset(isnan(classCurrent(pauseOnset,2))&classCurrent(pauseOnset,3)==2);
    pauseAftUnknBefUnkn = pauseOnset(isnan(classCurrent(pauseOnset,2))&isnan(classCurrent(pauseOnset,3)));
    
    %add to output structure
    windowGroup(19).edgeClassInfo = ...
        [windowGroup(19).edgeClassInfo; ...
        [pauseAftProtBefProt iWindow*ones(length(pauseAftProtBefProt),1) ...
        classCurrent(pauseAftProtBefProt,4:5)]];
    windowGroup(20).edgeClassInfo = ...
        [windowGroup(20).edgeClassInfo; ...
        [pauseAftProtBefRetr iWindow*ones(length(pauseAftProtBefRetr),1) ...
        classCurrent(pauseAftProtBefRetr,4:5)]];
    windowGroup(21).edgeClassInfo = ...
        [windowGroup(21).edgeClassInfo; ...
        [pauseAftProtBefUnkn iWindow*ones(length(pauseAftProtBefUnkn),1) ...
        classCurrent(pauseAftProtBefUnkn,4:5)]];
    windowGroup(22).edgeClassInfo = ...
        [windowGroup(22).edgeClassInfo; ...
        [pauseAftRetrBefProt iWindow*ones(length(pauseAftRetrBefProt),1) ...
        classCurrent(pauseAftRetrBefProt,4:5)]];
    windowGroup(23).edgeClassInfo = ...
        [windowGroup(23).edgeClassInfo; ...
        [pauseAftRetrBefRetr iWindow*ones(length(pauseAftRetrBefRetr),1) ...
        classCurrent(pauseAftRetrBefRetr,4:5)]];    
    windowGroup(24).edgeClassInfo = ...
        [windowGroup(24).edgeClassInfo; ...
        [pauseAftRetrBefUnkn iWindow*ones(length(pauseAftRetrBefUnkn),1) ...
        classCurrent(pauseAftRetrBefUnkn,4:5)]];    
    windowGroup(25).edgeClassInfo = ...
        [windowGroup(25).edgeClassInfo; ...
        [pauseAftUnknBefProt iWindow*ones(length(pauseAftUnknBefProt),1) ...
        classCurrent(pauseAftUnknBefProt,4:5)]];
    windowGroup(26).edgeClassInfo = ...
        [windowGroup(26).edgeClassInfo; ...
        [pauseAftUnknBefRetr iWindow*ones(length(pauseAftUnknBefRetr),1) ...
        classCurrent(pauseAftUnknBefRetr,4:5)]];    
    windowGroup(27).edgeClassInfo = ...
        [windowGroup(27).edgeClassInfo; ...
        [pauseAftUnknBefUnkn iWindow*ones(length(pauseAftUnknBefUnkn),1) ...
        classCurrent(pauseAftUnknBefUnkn,4:5)]];    
    
end

%% combination groups
windowGroup(9).edgeClassInfo = vertcat(windowGroup(1:8).edgeClassInfo);
windowGroup(18).edgeClassInfo = vertcat(windowGroup(10:17).edgeClassInfo);

%% ~~~ the end ~~~

 