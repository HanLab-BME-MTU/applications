function [roiIdx,mask] = createROI(MD,varargin)
% createROI - Displays MIPs of a volume and allows user to specify size of
% ROI.  Outputs the indices of the ROI boundaries, as well as a mask for
% uTrack.
%
% - Kevin Dean,2017.
%
% Required Inputs:
% MD               - Movie Data
%
% Optional Inputs:
% channel           - Which channel do you want to display?
% timeIdx           - What timepoint do you want to display?
%
% Outputs:
% roiIdx            - xmin, xmax, ymin, ymax, zmin, zmax
% mask              - binary mask with ones at ROI

ip = inputParser; ip.CaseSensitive = false;  ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('channel', 1, @isnumeric);
ip.addParameter('timeIdx', 1, @isnumeric);
ip.parse(MD, varargin{:});
p=ip.Results;

vol=double(MD.getChannel(p.channel).loadStack(p.timeIdx));
maxXY = squeeze(max(vol,[],3)); imshow(maxXY,[],'Border','tight'); xyIndices = ceil(getrect);

vol2 = vol(xyIndices(2):xyIndices(2)+xyIndices(4),xyIndices(1):xyIndices(1)+xyIndices(3),:);
maxYZ = squeeze(max(vol2,[],1)); imshow(maxYZ,[],'Border','tight'); zIndices = ceil(getrect);

roiIdx = nan(6,1);
roiIdx(1) = xyIndices(1);
roiIdx(2) = xyIndices(1)+xyIndices(3)-1;
roiIdx(3) = xyIndices(2);
roiIdx(4) = xyIndices(2)+xyIndices(4)-1;
roiIdx(5) = zIndices(1);
roiIdx(6) = zIndices(1)+zIndices(3)-1;

mask = zeros(size(vol,1),size(vol,2),size(vol,3));
mask(roiIdx(3):roiIdx(4),roiIdx(1):roiIdx(2),roiIdx(5):roiIdx(6)) = 1;

close all
