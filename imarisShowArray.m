function imarisApplication = imarisShowArray(array,imaApp)
%IMARISSHOWARRAY loads a 3D array into Imaris
%
% SYNOPSIS  imaApplication = imarisShowArray(array,imaApp)
%
% INPUT   array: Data to display. Usually 3D, but dimension can be lower.
%                It can be up to 5D - x,y,z,t,channel (this order).  
%         imaApp (opt): handle to the imaris application into which the array
%                       should be loaded. If no handle is given, Imaris
%                       will open a new session.
%
% OUTPUT  imarisApplication: Handle to the imaris Application
%
% REMARKS: (1) Imaris will close once all its handles are deleted. Even if
%              no output argument is assigned, a variable called
%              imarisApplication will be created in the workspace
%          (2) Imaris can not handle doubles. All values will be converted
%              to singles
%
% c: jonas, 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============
% TEST INPUT
%=============

% check dimensionality of input matrix
if nargin < 1 || isempty(array)
    error('you have to input at least the input array!')
end
% store size of array. Pad empty dimension with ones
sizeArray = ones(1,5);
% we have to permute the array, so that dimension 4 will be displayed in
% the time channel
array = permute(array,[1,2,3,5,4]);
sizeArray(1:ndims(array)) = size(array);
if length(sizeArray) > 5
    error('too many dimensions to display!')
end



% check whether we have to start a new imaris
if nargin < 2 || isempty(imaApp)
    imaApp = imarisStartNew;
else
    % make sure user did not accidentially close the window
    imaApp.mVisible = 1;
    % we do not want to assign a new imaris handle - take the old one!
    imaAppName = inputname(2);
end

% name of the input variable
varName = inputname(1);

%=============


%===============
% DISPLAY DATA
%===============

% create data set
imaDataSet = imaApp.mFactory.CreateDataSet;

% store matlab data
imaDataSet.SetData(single(array));

% set xyz-coordinates from 1:n
[imaDataSet.mExtendMinX,imaDataSet.mExtendMinY,imaDataSet.mExtendMinZ] = ...
    deal(1);
imaDataSet.mExtendMaxX = sizeArray(1);
imaDataSet.mExtendMaxY = sizeArray(2);
imaDataSet.mExtendMaxZ = sizeArray(3);

% do not set color table - is buggy

% adjust color range to min/max
minArray = min(array(:));
maxArray = max(array(:));

for ch = 0:size(array,4)-1
    imaDataSet.SetChannelRange(ch,minArray,maxArray);
    imaDataSet.SetChannelName(ch,num2str(ch+1));
end

% show in imaris
imaApp.mDataSet = imaDataSet;

% set view to projections
imaApp.mViewer = 'eViewerSection';

% name the imaris session
imaApp.mDataSet.SetParameter('Image', 'Name', varName);

%========================


%================
% ASSIGN OUTPUT
%================

if nargout > 0
    imarisApplication = imaApp;
end


%================