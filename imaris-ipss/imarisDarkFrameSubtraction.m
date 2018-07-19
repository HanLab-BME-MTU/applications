function correctedMovieOut = imarisDarkFrameSubtraction(numCorrectionFrames)
%IMARISDARKFRAMESUBTRACTION subtracts dark frames from a movie that is currently loaded in imaris
%
% SYNOPSIS imarisDarkFrameSubtraction(numCorrectionFrames)
%
% INPUT  numCorrectionFrames : (opt) [dark frames @ beginning of movie, d.f.@ end o.m.]
%
% OUPUT  correctedMovieOut : the movie, corrected
%
% c: jonas 04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===============
% TEST INPUT
%===============

% input test is performed later, we first have to talk to imaris to know
% what's going on

%===============



%================================
% CALL IMARIS & GET ALL THE INFO
%================================

vImarisApplication = actxserver('Imaris.Application');

vImarisDataSet = vImarisApplication.mDataSet;

% matlab does not (yet) allow us getting the whole movie from imaris, so we
% read it slice by slice. First, we need the size info, though
imarisMovieSize = [vImarisDataSet.mSizeX,...
        vImarisDataSet.mSizeY,...
        vImarisDataSet.mSizeZ,...
        vImarisDataSet.mSizeC,...
        vImarisDataSet.mSizeT];

% check whether there is a movie at all - else load a movie 
if ~all(imarisMovieSize)
    % load movie in imaris
    % FUTURE IMPROVEMENT - FILTERSPEC WITH IMARIS-FILETYPES
    [fileName,pathName] = uigetfile('*.*','Please select an image file');
    
    % handle user cancel
    if isequal(fileName,0)
        error('no movie loaded in IMARISDARKFRAMESUBTRACTION')
    end
        
    % load movie into imaris
    vImarisApplication.FileOpen([pathName fileName], 'reader=''All Formats''');
    
    % and get all the corresponding properties
    vImarisDataSet = vImarisApplication.mDataSet;
    imarisMovieSize = [vImarisDataSet.mSizeX,...
        vImarisDataSet.mSizeY,...
        vImarisDataSet.mSizeZ,...
        vImarisDataSet.mSizeC,...
        vImarisDataSet.mSizeT];

    if ~all(imarisMovieSize)
        error('Problem loading movie into Imaris')
    end

end

% check info on number of correctionFrames
if nargin == 0 | isempty(numCorrectionFrames)
    answer=inputdlg({'frames at the beginning of movie','frames at end of movie'},...
        'Please select the number of correction frames',1,{'5','5'});
    
    if isempty(answer)
        error('no correction frames selected in IMARISDARKFRAMESUBTRACTION')
    end
    
    % get number of correctionFrames
    numCorrectionFrames = str2double(answer)';
end

% take care of accidential negative or imaginary input
numCorrectionFrames = abs(numCorrectionFrames);

% test number of correction frames
if sum(numCorrectionFrames) == 0
    error('zero correction frames selected in IMARISDARKFRAMESUBTRACTION')
end
if sum(numCorrectionFrames) > imarisMovieSize(5)
    error('more correction frames than movie frames IMARISDARKFRAMESUBTRACTION')
end

% now read the whole movie from imaris
rawMovie = zeros(imarisMovieSize);
for t = 1:imarisMovieSize(5)
    for c = 1:imarisMovieSize(4)
        for z = 1:imarisMovieSize(3)
            rawMovie(:,:,z,c,t) = vImarisDataSet.GetDataSlice(z-1,c-1,t-1);
        end
    end
end

%========================


%========================
% CORRECT MOVIE
%========================

% function correctedMovie = darkFrameSubtraction(rawMovie,numCorrectionFrames)



% read correction image
correctionFrames = cat(5,rawMovie(:,:,:,:,1:numCorrectionFrames(1)),...
    rawMovie(:,:,:,:,end-numCorrectionFrames(2)+1:end));

% make one slice out of it (there could be a second color!)
% the squeeze takes care of the possibility that there is no second color
correctionSlice = mean(squeeze(mean(correctionFrames,5)),3);

% cut the beginning and the end off the raw movie
correctedMovie = rawMovie(:,:,:,:,numCorrectionFrames(1)+1:end-numCorrectionFrames);

% remember max intensity of movie
maxIntensity = max(rawMovie(:));

% we do not need the raw movie anymore
clear rawMovie

% remember new movie size
newMovieSize = size(correctedMovie);

% create a correctionMovie to subtract from correctedMovie
correctionMovie = repmat(correctionSlice,[1,1,newMovieSize(3:5)]);

% subtract correctionMovie
correctedMovie = correctedMovie - correctionMovie;

% set the min to zero and max to what it was  before
correctedMovie = (correctedMovie-min(correctedMovie(:)))*maxIntensity;

clear correctionMovie;
%----end function-----

imarisMovieSize = size(correctedMovie);

% Since putting back a movie into imaris resets all, we want to read as
% much image info as we can
imarisExtends = [vImarisDataSet.mExtendMaxX,...
        vImarisDataSet.mExtendMaxY,...
        vImarisDataSet.mExtendMaxZ];

for c = 1:imarisMovieSize(4)
    imarisColor(c)   = vImarisDataSet.GetColor(c-1);
end

% put corrected movie into imaris
vImarisDataSet.SetData(single(correctedMovie));

% set extents
vImarisDataSet.mExtendMaxX = imarisExtends(1);
vImarisDataSet.mExtendMaxY = imarisExtends(2);
vImarisDataSet.mExtendMaxZ = imarisExtends(3);

% set color
for c = 1:imarisMovieSize(4)
    vImarisDataSet.SetColor(imarisColor(c),c-1);
end

% assign output if necessary
if nargin == 1
    correctedMovieOut = correctedMovie;
end


