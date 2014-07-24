function imarisHandle = imarisSetMovie
%IMARISSETMOVIE loads a movie via Matlab into Imaris
%
% SYNOPSIS: imarisHandle = imarisSetMovie
%
% INPUT none
%
% OUTPUT imarisHandle: Handle to the current imaris session
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: Jonas Dorn
% DATE: 09-Feb-2006
%
%


% set maxSize to 100 MB
loadOpt.maxSize = 100000000;

% load movie
[movie, movieHeader, loadStruct] = ...
    cdLoadMovie('ask', [], loadOpt);

% launch Imaris
imarisHandle = imarisStartNew;

% create dataSet - set empty movie
imarisHandle.mDataSet.Create(3,movieHeader.numRows,movieHeader.numCols,...
    movieHeader.numZSlices,movieHeader.numWvs,movieHeader.numTimepoints);



% % assign size of movie
% imaDataSet = imarisHandle.mFactory.CreateDataSet;
% 
% imaDataSet.mSizeX = movieHeader.numCols;
% imaDataSet.mSizeY = movieHeader.numRows;
% imaDataSet.mSizeZ = movieHeader.numZSlices;
% imaDataSet.mSizeC = movieHeader.numWvs;
% imaDataSet.mSizeT = movieHeader.numTimepoints;

% assign pixelsize (when will they allow that directly??)
% set xyz-coordinates 0:n*pixelSize
imarisHandle.mDataSet.mExtendMinX = 0;
imarisHandle.mDataSet.mExtendMinY = 0;
imarisHandle.mDataSet.mExtendMinZ = 0;
imarisHandle.mDataSet.mExtendMaxX = imarisHandle.mDataSet.mSizeX * movieHeader.pixelX;
imarisHandle.mDataSet.mExtendMaxY = imarisHandle.mDataSet.mSizeY * movieHeader.pixelY;
imarisHandle.mDataSet.mExtendMaxZ = imarisHandle.mDataSet.mSizeZ * movieHeader.pixelZ;

% now loop through movie and set volumes into DataSet
done = 0;
while ~done

    for t = 1:length(loadStruct.loadedFrames)
        currentTime = loadStruct.loadedFrames(t);
        for c = 1:1 % make color work later
            % put volume            
            imarisHandle.mDataSet.SetDataVolume(single(movie(:,:,:,c,t)),c-1,currentTime-1);
        end % color
    end % time

    % load move movie
    if isempty(loadStruct.frames2load)
        done = 1;
    else
        movie = [];
        [movie, movieHeader, loadStruct] = ...
            cdLoadMovie(loadStruct.movieType,[],loadStruct);
    end
end

% finish

% name the imaris session
imarisHandle.mDataSet.SetParameter('Image', 'Name', loadStruct.movieName);