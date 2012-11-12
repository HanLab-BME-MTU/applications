function [ fcorr] = correctDrift(features,drift,varargin)
%CORRECTDRIFT corrects localizations by given drift array
%
%   Microscope stage drift is a plague, so one has to correct for it. By
%   using fiduciary markers one can track the stage movement over time and
%   correct for these deviations. CORRECTDRIFT takes two arguments: 
%    * the detected localization as obtained from analyzeLMdata
%    * the calculated drift from getDrift
%   The output are the drift-corrected localizations in fcorr whose
%   internal structure is the same as the input cell array features.
%
%   INPUT:
%       features  ->  cell array as obtained from analyzeLMdata
%          drift  ->  drift in xy between two frames
%
%   OUTPUT:
%       fcorr  ->  localizations corrected for drift
%
%  US, 2012/10/25
%

ip=inputParser;

ip.CaseSensitive=true;
ip.StructExpand=true;

ip.addRequired('features',@iscell);
ip.addRequired('drift',@isnumeric);

ip.parse(features,drift,varargin{:});

F=ip.Results.features;
drift=ip.Results.drift;

nFrames=numel(F);

% make sure one has drift correction for every time point 
if numel(drift(:,1)) ~= nFrames
    msg='dimension mismatch between features and drift array';
    error(msg);
end

% cumulative drift from reference point
drift=cumsum(drift);

fcorr=F;
for k=1:nFrames
    
    if ~isempty(fcorr{k})
        pos=[fcorr{k}.x' fcorr{k}.y'];
        nLoc=numel(pos(:,1));
    
        for i=1:nLoc
            pos(i,:)=pos(i,:)-drift(k,:);
        end
    
        fcorr{k}.x=pos(:,1);
        fcorr{k}.y=pos(:,2);
    end
        
end



end

