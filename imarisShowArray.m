function imarisApplication = imarisShowArray(array,imaApp,pixelSize,equalizeTC,assigninBase)
%IMARISSHOWARRAY loads a 3D array into Imaris
%
% SYNOPSIS  imaApplication = imarisShowArray(array,imaApp,pixelSize,equalizeTC)
%
% INPUT   array: Data to display. Usually 3D, but dimension can be lower.
%                It can be up to 5D - x,y,z,t,channel (this order).  
%         imaApp (opt): handle to the imaris application into which the array
%                       should be loaded. If no handle is given, Imaris
%                       will open a new session.
%         pixelSize (opt) : [x,y,z] of pixel lengths {[1,1,1]}
%         equalizeTC (opt): [eqT, eqC] - switches indicating whether to
%                           equalize values along the 4th and/or 5th
%                           dimension. {[0,0]}
%         assigninBase (opt): whether to assign an ImarisHandle in the base
%                               workspace {1}
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

def_pixelSize = [1,1,1];
def_equalize = [0,0];
def_assigninBase = true;

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

if nargin < 5 || isempty(assigninBase)
    assigninBase = def_assigninBase;
end

%Setup colors for channels
if sizeArray(4) <= 3
    chanCols = [ 1 0 0 ; 0 1 0 ; 0 0 1];
else
    chanCols = jet(sizeArray(4));
end

% check whether we have to start a new imaris
if nargin < 2 || isempty(imaApp)
    try
        imaApp = imarisStartNew(assigninBase);
        is7pt5orLater = false;
    catch err        
        if strcmp(err.identifier,'MATLAB:COM:InvalidProgid')
            %If we need the newer version, use the IceImarisConnector
            is7pt5orLater = true;
            iceConn = IceImarisConnector;            
            iceConn.startImaris;
            imaApp = iceConn.mImarisApplication;
        else
            throw(err)
        end
    end
        
else
    % make sure user did not accidentially close the window
    imaApp.mVisible = 1;
end

% test which pixelsize to use
if nargin < 3 || isempty(pixelSize)
    pixelSize = def_pixelSize;
else
    pixelSize = returnRightVector(pixelSize, 3);
end

if nargin < 4 || isempty(equalizeTC)
    equalizeTC = def_equalize;
else
    if length(equalizeTC) == 1
        equalizeTC = [equalizeTC, equalizeTC];
    end   
end

% name of the input variable
varName = inputname(1);
%=============

%======================
% EQUALIZE IF SELECTED
%======================

switch sum(find(equalizeTC == 1))
    % loop to avoid memory problems
    case 1 % only t
        for t=1:sizeArray(4)
            frame = array(:,:,:,t,:);
            frameVec = frame(:);
            minFrame = nanmin(frameVec);
            maxFrame = nanmax(frameVec);
            array(:,:,:,t,:) = ...
                (frame - minFrame) / (maxFrame-minFrame);
        end
            
    case 2 % only c
        for c=1:sizeArray(5)
            frame = array(:,:,:,:,c);
            frameVec = frame(:);
            minFrame = nanmin(frameVec);
            maxFrame = nanmax(frameVec);
            array(:,:,:,:,c) = ...
                (frame - minFrame) / (maxFrame-minFrame);
        end
    case 3 % t and c
        for c=1:sizeArray(5)
            for t=1:sizeArray(4)
                frame = array(:,:,:,t,c);
            frameVec = frame(:);
            minFrame = nanmin(frameVec);
            maxFrame = nanmax(frameVec);
            array(:,:,:,t,c) = ...
                (frame - minFrame) / (maxFrame-minFrame);
            end
        end
        
    otherwise % don't equalize
        
end

% get rid of NaNs
delta = 0;
for c = 1:sizeArray(5)
    for t = 1:sizeArray(4)
       frame = array(:,:,:,t,c);
            frameVec = frame(:);
            nanList = (isnan(frameVec));
            if any(nanList)
                if delta == 0
                    delta = nanmin(diff(unique(frameVec)));
                end
                frame(nanList) = nanmin(frameVec)-delta;
                array(:,:,:,t,c) = frame;
            end
    end
end

%===============
% DISPLAY DATA
%===============

if ~is7pt5orLater

    % create data set
    imaDataSet = imaApp.mFactory.CreateDataSet;

    % store matlab data
    imaDataSet.SetData(single(array));

    % set xyz-coordinates from 1:n or pixelSize:n*pixelSize
    imaDataSet.mExtendMinX = pixelSize(1);
    imaDataSet.mExtendMinY = pixelSize(2);
    imaDataSet.mExtendMinZ = pixelSize(3);
    imaDataSet.mExtendMaxX = sizeArray(1) * pixelSize(1);
    imaDataSet.mExtendMaxY = sizeArray(2) * pixelSize(2);
    imaDataSet.mExtendMaxZ = sizeArray(3) * pixelSize(3);

    % do not set color table - is buggy

    % adjust color range to min/max
    minArray = nanmin(array(:));
    maxArray = nanmax(array(:));

    for ch = 0:size(array,4)-1
        imaDataSet.SetChannelColor(ch,chanCols(ch+1,1),chanCols(ch+1,2),chanCols(ch+1,3),0);
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
else
    %For newer imaris versions using the IceConnector, the commands are
    %(annoyingly) slightly different
    
    % create data set
    imaDataSet = imaApp.GetFactory.CreateDataSet;
    classDataSet=Imaris.tType.eTypeFloat;
    imaDataSet.Create(classDataSet ,sizeArray(1),sizeArray(2),sizeArray(3),sizeArray(4),sizeArray(5)); 
    imaApp.SetDataSet(imaDataSet);

    % store matlab data
    for ch = 1:size(array,4)
        for  t = 1:size(array,5)            
            iceConn.setDataVolume(single(array(:,:,:,ch,t)),ch-1,t-1);
            pause(.5);%For whatever reason if you add the next channel too quickly there are some weird glitches. Stupid temporary fix.
        end
    end    

    % set xyz-coordinates from 1:n or pixelSize:n*pixelSize
    imaDataSet.SetExtendMinX(pixelSize(1));
    imaDataSet.SetExtendMinY(pixelSize(2));
    imaDataSet.SetExtendMinZ(pixelSize(3));
    imaDataSet.SetExtendMaxX(sizeArray(1) * pixelSize(1));
    imaDataSet.SetExtendMaxY(sizeArray(2) * pixelSize(2));
    imaDataSet.SetExtendMaxZ(sizeArray(3) * pixelSize(3));
    
    % do not set color table - is buggy

    % adjust color range to min/max
    minArray = nanmin(array(:));
    maxArray = nanmax(array(:));

    for ch = 0:size(array,4)-1
        imaDataSet.SetChannelColorRGBA(ch,iceConn.mapRgbaVectorToScalar([chanCols(ch+1,:),0]));
        imaDataSet.SetChannelRange(ch,minArray,maxArray);
        imaDataSet.SetChannelName(ch,num2str(ch+1));
    end
        
    % name the imaris session
    imaDataSet.SetParameter('Image', 'Name', varName);

    %Adjust the view
    imaApp.GetSurpassCamera.Fit;    
    
end

%================
% ASSIGN OUTPUT
%================

if nargout > 0
    imarisApplication = imaApp;
end


%================