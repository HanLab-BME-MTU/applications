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

% assign size of movie
imaDataSet = imarisHandle.mFactory.CreateDataSet;

imaDataSet.mSizeX = movieHeader.numCols;
imaDataSet.mSizeY = movieHeader.numRows;
imaDataSet.mSizeZ = movieHeader.numZSlices;
imaDataSet.mSizeC = movieHeader.numWvl;
imaDataSet.mSizeT = movieHeader.numTimepoints;

% assign pixelsize (when will they allow that directly??)
% set xyz-coordinates 0:n*pixelSize
imaDataSet.mExtendMinX = 0;
imaDataSet.mExtendMinY = 0;
imaDataSet.mExtendMinZ = 0;
imaDataSet.mExtendMaxX = imaDataSet.mSizeX * movieHeader.pixelX;
imaDataSet.mExtendMaxY = imaDataSet.mSizeY * movieHeader.pixelY;
imaDataSet.mExtendMaxZ = imaDataSet.mSizeZ * movieHeader.pixelZ;

% now loop through movie and set volumes into DataSet
done = 0;
while ~done

    for t = 1:loadStruct.loadedFrames
        for c = 1:1 % make color work later
            % put volume
            imaDataSet.SetVolume(single(movie(:,:,:,c,t)),c-1,t-1);
        end % color
    end % time

    % load move movie
    if ~isempty(loadStruct.frames2load)
        done = 1;
    else
        movie = [];
        [movie, movieHeader, loadStruct] = ...
            cdLoadMovie(loadStruct);
    end
end

% finish

% show in imaris
imaApp.mDataSet = imaDataSet;

% name the imaris session
imaApp.mDataSet.SetParameter('Image', 'Name', loadStruct.movieName);