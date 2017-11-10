function trackSEL = getSpeckleTrackSEL(MD,varargin)
% GETSPECKLETRACKSEL outputs start times, end times and lifetimes of speckles
%
%SYNOPSIS trackSEL = getSpeckleTrackSEL(MD);
%
%INPUT  
%       MD - a MovieData object
%
%       iChan - (optional) the channel containing the tracks. If not input
%       the first channel containing valid speckle tracks will be used.
%
%OUTPUT 
%      trackSEL          : An array with 3 columns and number of rows equal
%                           to number of (compound) tracks. 1st column 
%                           indicates track start times, 2nd column
%                           indicates track end times and 3rd column
%                           indicates track lifetimes.
%
% Note: should be consistent with getTrackSEL (for Khuloud tracks
% structure)

% Sebastien Besson, Jan 2012

% Input check
ip = inputParser;
ip.addRequired('MD',@(x) isa(x,'MovieData'));
ip.addOptional('iChan',[],@isscalar);
ip.parse(MD,varargin{:});
iChan=ip.Results.iChan;

% Retrieve index of speckle tracking process
iProc = MD.getProcessIndex('SpeckleTrackingProcess',1,false);
assert(~isempty(iProc),'No speckle tracking process found');

% Check there is an output
specTrackProc=MD.processes_{iProc};
assert(specTrackProc.success_ ,'Please run speckle tracking first');

% Load tracking output
if isempty(iChan)
    assert(any(specTrackProc.checkChannelOutput),'No valid output found');
    iChan = find(specTrackProc.checkChannelOutput,1);
else
    assert(specTrackProc.checkChannelOutput(iChan),...
        'Input channel does not have valid speckle tracks');
end

% Load the MPM and keep x-components only
MPM=specTrackProc.loadChannelOutput(iChan,'output','MPM');
MPM = MPM(:,1:2:end)';
% Add zeros a time nFrames+1
MPM(end+1,:)=0;

% Retrieve linear index of non-tracks elemetns
nontrackTimes=find([1; ~MPM(:); 1]);
% Compute lifetime of block betwen two consecutive non-track elements
lifetimes = diff(nontrackTimes)-1;

% Determine number of tracks and lifetimes
nTracks=sum(lifetimes~=0);
trackSEL = zeros(nTracks,3);
trackSEL(:,3) = lifetimes(lifetimes~=0);

% Compute start and end of tracks
% Linear index is converted back into time index
% Should substract 1 as we added an extra 1 at the beginning of the array
% Should add one as we are interest in the start of the track (i.e. the
% next element after the nontrackTimes)
trackSEL(:,1)  =  mod(nontrackTimes(lifetimes~=0),size(MPM,1));
trackSEL(:,2)  =  trackSEL(:,1)+trackSEL(:,3) -1;

assert(min(trackSEL(:,1))>=1 && max(trackSEL(:,1))<=MD.nFrames_);
assert(max(trackSEL(:,3))<=MD.nFrames_);
end