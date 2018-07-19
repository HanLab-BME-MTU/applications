function imarisBleachingCorrection(imarisAppOrID)
%IMARISBLEACHINGCORRECTION corrects a movie for bleaching with an exponential model for the intensity decay
%
% SYNOPSIS imarisBleachingCorrection(imarisAppOrID)
%
% INPUT    imarisAppOrID : handle to Imaris Application or Application ID
%
% The program will fit a model of the form a*exp(-t/b)+c to the intensity
% of the movie data of Imaris. This is based on the assumption that the
% amount of fluorophores in the image is constant or growing much more
% slowly than the decay (other cells moving into the image will adversely
% affect the calculation!)
% 
% c: jonas 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% TEST INPUT
%============

try
    imaApp = imarisCheckHandle(imarisAppOrID);
catch
    error('please specify a valid imaris handle or id')
end

%============

%===================
% READ IMARIS DATA
%===================

% handle to the data
imaDataSet = imaApp.mDataSet;

% data size
dataSize = zeros(1,5);
dataSize(1) = imaDataSet.mSizeX;
dataSize(2) = imaDataSet.mSizeY;
dataSize(3) = imaDataSet.mSizeZ;
dataSize(4) = imaDataSet.mSizeC;
dataSize(5) = imaDataSet.mSizeT;

% read movie
movie = zeros(dataSize);

for c=1:dataSize(4)
    for t = 1:dataSize(5)
        for z = 1:dataSize(3)
            % Remember -1
            movie(:,:,z,c,t) = imaDataSet.GetDataSlice(z-1,c-1,t-1);
        end
    end
end

%==================

%=====================
% SELECT MODE/CHANNEL
%=====================

% later: inputGUI mit choice of channel, choice of slicewise or framewise

selectedChannel = 1;
mode = 'frame';

%====================

%====================
% ESTIMATE BLEACHING
%====================

% current model: A*exp(-t/B)+C

% read intensities
switch mode
    case 'frame'
        intensity = squeeze(...
            sum(sum(sum(movie(:,:,:,selectedChannel,:),1),2),3));
        % start x at zero: in the beginning, there was no bleaching
        xData = 0:dataSize(5)-1;
        
    case 'slice'
        
        error('slice mode not implemented yet')
        
        intensity = ...
            sum(sum(movie(:,:,:,selectedChannel,:),1),2);
        xLength = dataSize(5)*dataSize(3);
        intensity = reshape(intensity,[xLength,1]);
        
        % start x at zero: in the beginning, there was no bleaching
        xData = 0:xLength-1;
end

% fit with robust tools
expFit = robustExponentialFit(xData,intensity,1);

% show data
figure('Name','Exponential Fit');
yExp = expFit(1)*exp(-xData/expFit(2))+expFit(3);
plot(xData,intensity,'.k',xData,yExp,'r');


%====================

        
%=====================
% CORRECT
%=====================

% 1) subtract background
% 2) divide by estimate
% 3) make channel [0,1]

switch mode
    case 'frame'
        background = expFit(3)/prod(dataSize(1:3));
        bleachingEstimate = yExp - expFit(3);
        
        % in order to save memory, we divide the frames individually
        for t = 1:dataSize(4)
            movie(:,:,:,t,selectedChannel) = ...
                movie(:,:,:,selectedChannel,t)./bleachingEstimate(t);
        end
        
        % to make [0,1]: minus minimum, divide by maximum
        movie(:,:,:,selectedChannel,:) = movie(:,:,:,selectedChannel,:) -...
            min(reshape(movie(:,:,:,selectedChannel,:),1,[]));
        movie(:,:,:,selectedChannel,:) = movie(:,:,:,selectedChannel,:) ./...
            max(reshape(movie(:,:,:,selectedChannel,:),1,[]));
        
    case 'slice'
        error('not implemented yet')
end

%=====================

%=====================
% WRITE BACK
%=====================

% set undo
imaApp.DataSetPushUndo('BleachingCorrection')
% write back movie
imaDataSet.SetData(single(movie));
        
