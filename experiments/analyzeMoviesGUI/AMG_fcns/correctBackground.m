function [newDataMovie,newDataMovieName,newOldCorrectionInfo] = correctBackground(dataMovieName,correctionInfo,oldCorrectionInfo,fid1,fid2)
%correctBackground makes a subtraction of the background from the dataMovie according to correctionInfo. 

%load dataMovie. We're already in the correct directory, and the user can not
%change it due to the waitfor command

%check input: fids for logging to file
if nargin < 4
    fid1 = [];
    fid2 = [];
elseif nargin < 5 
    if ~isempty(fid1)
        fid2 = 0;
    else
        fid2 = [];
    end
end

%make sure we're loading a r3d-movie
if strcmp(dataMovieName(end),'c')
    %this could either be a cropped movie or a corrected movie
    if strcmp(dataMovieName(end-8:end),'_crop.r3c')
        %the movie has been cropped, ergo it has been corrected already (or it
        %is old and hasn't, but that's not interesting here, anyway)
        
        %the calling programs are waiting for a movie and a name - let's hand it
        %to them
        newDataMovie = readmat(dataMovieName);
        newDataMovieName = dataMovieName;
        newOldCorrectionInfo = oldCorrectionInfo;
        return
        
    elseif strcmp(dataMovieName(end-8:end),'_corr.r3c')
        %the movie has been corrected before
        if isequalwithequalnans(correctionInfo,oldCorrectionInfo)
            %the correctionInfo has not changed, so we can use the same movie
            newDataMovie = readmat(dataMovieName);
            newDataMovieName = dataMovieName;
            newOldCorrectionInfo = oldCorrectionInfo;
            return
        else
            %the correctionInfo has changed - not good
            errordlg('you can''t correct a cropped movie')
            return
        end
    else
        %this is some strange moviename - not supported
        error('check dataMovieName!')
    end
end

if fid1
    fprintf(fid1,[nowString,' movie  =  r3dread(%s);\n'],dataMovieName);
    fprintf(fid2,[nowString,' load unfiltered movie: (%s)\n'],dataMovieName);
end

%read original movie
dataMovie = r3dread(dataMovieName);

%assign newOldCorrectionInfo
newOldCorrectionInfo = correctionInfo;

%make sure there is a correctionInfo
if isempty(correctionInfo)
    %there is no correctionInfo: return r3d-movie
    newDataMovie = dataMovie;
    newDataMovieName = dataMovieName;
    newOldCorrectionInfo = oldCorrectionInfo;
    return
end


%switch on correction
switch isempty(correctionInfo.correctFrames)+2*isempty(correctionInfo.header)
    %good outcomes: 0, 1; otherwise the header is missing ->error
    case 1 %correction with darkImage-movie(s)
        oldDir = pwd;
        corrImg = [];
        for i = 1:size(correctionInfo.correctFiles,1)
            %move to corrMovieDir via bioDataMainDir
            cdBiodata(0);
            cd(correctionInfo.correctFiles{i,1});
            
            if fid1
                fprintf(fid1,[nowString,' corrMov  =  r3dread(%s\%s);\n'],correctionInfo.correctFiles{i,1},correctionInfo.correctFiles{i,2});
                fprintf(fid2,[nowString,' load correction movie: (%s\%s)\n'],dataMovieName);
            end
            
            %read correctionMovie
            corrMovieName = correctionInfo.correctFiles{i,2};
            if strcmp(corrMovieName(end),'d');
                corrMov = r3dread(corrMovieName);
            else
                corrMov = readmat(corrMovieName);
            end
            %average over timepoints
            corrImgTmp = squeeze(mean(corrMov,5));
            %delete correctionMove
            clear corrMov;
            %store corrImg (max 2 movies, so no problem with mean)
            corrImg = squeeze(mean(cat(4,corrImg,corrImgTmp),4));
            
        end
        %store only one slice
        corrImg = mean(corrImg,3);
        
        %return to movieDir
        cd(oldDir);
        %get movieSize
        movieSize = size(dataMovie);
        
    case 0 %correction with firs/last X frames
        %cut first X and last Y frames (this command also works if one entry of cF is 0)
        
        if fid1
            fprintf(fid1,[nowString,' corrFrames  =  dataMovie(:,:,:,:,[1:%i,end-%i+1:end]);\n'],correctionInfo.correctFrames(1),correctionInfo.correctFrames(2));
            fprintf(fid2,[nowString,' extract correction frames: first %s, last %s)\n'],correctionInfo.correctFrames(1),correctionInfo.correctFrames(2));
        end
        
        corrFrames = dataMovie(:,:,:,:,[1:correctionInfo.correctFrames(1),end-correctionInfo.correctFrames(2)+1:end]);
        %shorten dataMovie
        dataMovie = dataMovie(:,:,:,:,correctionInfo.correctFrames(1)+1:end-correctionInfo.correctFrames(2));
        %calculate corrImg
        corrImg = squeeze(mean(corrFrames,5));
        
        %store only one slice
        corrImg = mean(corrImg,3);
        
        %correct header.time. Since it is a 1 by numzSlices*numTimepoints
        %vector, we have to calculate from where to where to take it
        movieSize = size(dataMovie);
        firstNum = (correctionInfo.correctFrames(1))*movieSize(3)+1; 
        lastNum = (movieSize(5)+correctionInfo.correctFrames(1))*movieSize(3); %ms5+cF1+cF2-cF2
        hTime = correctionInfo.header.Time(firstNum:lastNum);
        
    otherwise
        error('movie header is missing')
end



%store correctionData (don't store the header twice!)
correctionData.info.correctFrames = correctionInfo.correctFrames;
correctionData.info.correctFiles  = correctionInfo.correctFiles;
correctionData.image = corrImg;

if fid1
    fprintf(fid1,[nowString,' save(''correctionData'',''correctionData'');\n']);
    fprintf(fid2,[nowString,' store correctionData\n']);
end

save('correctionData','correctionData');

%get movieLength
numTimepoints = movieSize(5);
numZSlices = movieSize(3);

%make a correctionMovie
correctionMovie = repmat(corrImg,[1,1,numZSlices,1,numTimepoints]);

%create new DataMovie
newDataMovie = dataMovie-correctionMovie;

clear correctionMovie;

%make sure there are no negative values in the movie
newDataMovie = newDataMovie - min(newDataMovie(:));

% %normalize Movie to camera data range (12-bit-data!)
% newDataMovie = newDataMovie/MOVIERANGE;

%create new fileName
newDataMovieName = [dataMovieName(1:end-4),'_corr.r3c'];

%make sure there is not already a corrected movie
if exist(newDataMovieName)
    if fid1
        fprintf(fid1,[nowString,' delete(%s);\n'],newDataMovieName);
        fprintf(fid2,[nowString,' remove old corrected movie %s\n'],newDataMovieName);
    end
    delete(newDataMovieName);
end

if fid1
    fprintf(fid1,[nowString,' writemat(%s, newDataMovie);\n'],newDataMovieName);
    fprintf(fid2,[nowString,' save movie %s\n'],newDataMovieName);
end

%write newDataMovie to file
writemat(newDataMovieName, newDataMovie);

%read header
r3dMovieHeader = correctionInfo.header;

%correct time, movieLength
r3dMovieHeader.numTimepoints = numTimepoints;
if exist('hTime','var')
    r3dMovieHeader.Time = hTime;
end
r3dMovieHeader.cropInfo = []; %for compatibility
r3dMovieHeader.correctInfo = correctionData.info;

if fid1
    fprintf(fid1,[nowString,' save(''r3dMovieHeader'',''r3dMovieHeader'');\n']);
    fprintf(fid2,[nowString,' save updated movie header\n']);
end

%save movieHeader
save('r3dMovieHeader','r3dMovieHeader');

