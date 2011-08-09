function angles =plusTipPlotTrackAngles(projData,saveDir)
% plusTipPlotTrackAngles read the project tracks, get the angles and plot them
%
% Input:
% 
%   projData  : structure containing field .anDir, which gives the full 
%               path to the roi_x directory 
%

% Sebastien Besson, August 2011

% Input check
ip =inputParser;
ip.addRequired('projData',@isstruct)
ip.addRequired('saveDir',@ischar);
ip.parse(projData,saveDir);

% Load tracking result
trackDir = [projData.anDir filesep 'track'];
try
    files = dir([trackDir filesep '*.mat']);
    s=load([trackDir filesep files(1).name],'tracksFinal');
catch ME
    error(['--plusTipPostTracking:' ME.message]);
end

% Get all track angles and concatenate them
angles = arrayfun(@getAngles,s.tracksFinal,'UniformOutput',false);
angles =horzcat(angles{:});
angles=angles(~isnan(angles));


saveFig =figure;
rose(angles);
print(saveFig,'-dtiff',[saveDir filesep 'angles_histogram.tif'])
close(saveFig)

end

function angles = getAngles(track)
x=track.tracksCoordAmpCG(1:8:end);
y=track.tracksCoordAmpCG(2:8:end);
angles=atan2(y(2:end)-y(1:end-1),x(2:end)-x(1:end-1));
end