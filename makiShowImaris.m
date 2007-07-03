function makiShowImaris(dataStruct,select)
%MAKISHOWIMARIS shows mammalina kinetochore data in Imaris
%
% SYNOPSIS: imarisHandle = makiShowImaris(dataStruct,select)
%
% INPUT dataStruct: data structure as described in makiMakeDataStruct
%                   if empty, guiLoad
%		select: selection switch.
%               1: show initCoord (default)
%
% OUTPUT ---
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 29-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input
if nargin < 1
    dataStruct = [];
end
if nargin < 2 || isempty(select)
    select = 1;
end
if isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

% turn off property reader warning
warningState = warning;
warning off IMARISIMREAD:NOPROPERTYREADER

% reduce amount of typing
dataProperties = dataStruct.dataProperties;
pixelSize = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];
% check for croppping
if isempty(dataProperties.crop)
    dataProperties.crop = zeros(2,3);
end
crop = dataProperties.crop(:,1:3);
isCrop = any(crop,1);
crop(1,~isCrop) = 1;
crop(2,~isCrop) = dataProperties.movieSize(find(~isCrop)); %#ok<FNDSB>


% loop switches

for sw = select(:)'
    
    % start imaris
            imarisApplication = imarisStartNew;
            % load raw movie into imaris. We could do filtered movie, but
            % for this, we would have to load frame by frame and do all the
            % image properties stuff
            imarisApplication.FileOpen(...
                fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName),...
                'reader=''DeltaVision''');
            
            % check image properties: image should begin at -0.5 pix.
            % It would be nice to be able to set the pixelSize to -0.5.
            % Instead, we have to read the mins and calculate an offset
            zeroOffsetX = imarisApplication.mDataSet.mExtendMinX + 0.5*pixelSize(1);
            zeroOffsetY = imarisApplication.mDataSet.mExtendMinY + 0.5*pixelSize(2);
            zeroOffsetZ = imarisApplication.mDataSet.mExtendMinZ + 0.5*pixelSize(3);
            zeroOffset = [zeroOffsetX zeroOffsetY zeroOffsetZ];
    
    
    switch sw
        case 1
            
            
            
            % plot from initCoord

            initCoord = dataStruct.initCoord;

            % make spots object: X,Y,Z,T,r

            nTimepoints = dataProperties.movieSize(end);
            nSpots = zeros(nTimepoints,1);
            for t=1:nTimepoints
                nSpots(t) = size(initCoord{t},1);
            end
            spots = zeros(sum(nSpots),5);
            spots(:,5) = pixelSize(1)*2; % radius in micron
            goodTimes = find(nSpots);
            nSpotSum = [0;cumsum(nSpots)];
            for t = goodTimes'
                % calculate positions in microns. Subtract one voxel from
                % the coords: Imaris starts counting at 0!
                
                spots(nSpotSum(t)+1:nSpotSum(t+1),1:4) = ...
                    [(cat(1,initCoord{t}(:,[2,1,3]))-1 + ...
                    repmat(crop(1,[2,1,3])-1,nSpots(t),1)).*...
                    repmat(pixelSize,nSpots(t),1) + ...
                    repmat(zeroOffset,nSpots(t),1),...
                    (t-1)*ones(nSpots(t),1)];
            end

            
             
            


            % make top-level surpass scene
            imaSurpassScene = imarisApplication.mFactory.CreateDataContainer();

            % fill surpass scene with light and frame and volume
            imaLight = imarisApplication.mFactory.CreateLightSource();
            imaSurpassScene.AddChild(imaLight);
            imaFrame = imarisApplication.mFactory.CreateFrame();
            imaSurpassScene.AddChild(imaFrame);
            imaVolume = imarisApplication.mFactory.CreateVolume();
            imaSurpassScene.AddChild(imaVolume);

            % add surpass scene and set view
            imarisApplication.mSurpassScene = imaSurpassScene;
            imarisApplication.mViewer = 'eViewerSurpass';

            % create spots object
            imaSpots = imarisApplication.mFactory.CreateSpots;

            % set coords
            imaSpots.Set(single(spots(:,1:3)),single(spots(:,4)),single(spots(:,5)));

            % add to scene
            imaSurpassScene.AddChild(imaSpots);

        otherwise
            fprintf(1,'selection %i not implemented yet',sw)
    end
end

% turn warnings back on
warning(warningState)