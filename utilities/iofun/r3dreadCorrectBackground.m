function [backgroundCorrectionImage,minStart,maxEnd,r3dMovieHeader] = r3dreadCorrectBackground(filename,numRow,numCol,numZ,numWvl,numTimes,correctionDataFromFile,correctionDataFromVar,MOVIERANGE,r3dMovieHeader)
%CORRECTBACKGROUND calculates the backgroundcorrection image for movies loaded with r3dread
%
% INPUT filename: moviefilename
%       num...  : movieSize
%       correctionData... : information about how to correct
%       MOVIERANGE : constant for calling r3dread
%
% OUTPUT image of size [numRow,numCol,numZ,numWvl] describing the
%        background.
%        minStart: frame from which on the movie should be read
%        maxEnd  : last interesting frame
%
% The program will also save correctionData and r3dMovieHeader to disk, if
% a new correctionData has been calculated.


% First: find out whether we already have a correctionImage
% if there is a correctionFile and the info is the same as in the var, we
% alreay have an image. If the correctionInfo of the var is 'file', we can
% go ahead and calculate the image already. Else we call r3dread for the
% first and last n frames and calculate the BG image - we want to be able
% to save RAM!!


    isThereVar = ~isempty(correctionDataFromVar);
    isThereFile = ~isempty(correctionDataFromFile);
    
    if isThereVar % change var
        info = struct('correctFiles',correctionDataFromVar.correctFrames,...
            'correctFrames',correctionDataFromVar.correctFrames);
        correctionDataFromVar.info = info;
        correctionDataFromVar = rmfield(correctionDataFromVar,{'correctFiles','correctFrames'});
    end
    
    switch isThereFile + 2*isThereVar
        
        case 1 % file only, no var : load corrImg
            
            backgroundCorrectionImage = correctionDataFromFile.image;
            calcCorrection = 0;
            correctionData = correctionDataFromFile;
            
        case 2 % var only : set calcCorrection
            
            if isempty(correctionDataFromVar.info.correctFrames)
                calcCorrection = 2; % calcFromFile
            else
                calcCorrection = 1; % calcFromFrames
            end
            correctionData = correctionDataFromVar;
            
        case 3 % both. Check whether we have to do calcCorrection
            
            % if there info is the same, we are up to date
            
            if isequalwithequalnans(correctionDataFromVar.info,correctionDataFromFile.info)
                backgroundCorrectionImage = correctionDataFromFile.image;
                calcCorrection = 0;
                correctionData = correctionDataFromFile;
            else % this should be pretty rare...
                if isempty(correctionDataFromVar.info.correctFrames)
                    calcCorrection = 2; % calcFromFile
                else
                    calcCorrection = 1; % calcFromFrames
                end
                correctionData = correctionDataFromVar;
            end
            
        otherwise % should not happen, because correctBG should be 0
            backgroundCorrectionImage = zeros(numRow,numCol,numZ,numWvl);
            correctionData = [];
    end
    
    % now calculate backgroundCorrection
    switch calcCorrection
        
        case 0 % make sure we have just one slice!
            
            if size(backgroundCorrectionImage,3)~=1
                backgroundCorrectionSlice = mean(backgroundCorrectionImage,3);
            else
                backgroundCorrectionSlice = backgroundCorrectionImage;
            end
            
            backgroundCorrectionImage = repmat(backgroundCorrectionSlice,[1,1,numZ,1]);
            
            
        case 1 % calculate from frames
            
            %load first and last frames via r3dread
            firstCorrFrames = r3dread(filename,1,correctionData.info.correctFrames(1),MOVIERANGE,0);
            lastCorrFrames  = r3dread(...
                filename,numTimes-correctionData.info.correctFrames(2)+1,correctionData.info.correctFrames(2),MOVIERANGE,0);
            corrFrames = cat(5,firstCorrFrames,lastCorrFrames);
            
            % make correctionSlice
            corrFrameSingle = mean(corrFrames,5);
            backgroundCorrectionSlice = mean(corrFrameSingle,3);
            
            % this is compatible with several wavelengths
            backgroundCorrectionImage = repmat(backgroundCorrectionSlice,[1,1,numZ,1]);
            
            
            
        case 2 % calculate from file
            
            oldDir = pwd;
            corrImg = [];
            for i = 1:size(correctionData.info.correctFiles,1)
                %move to corrMovieDir via bioDataMainDir
                cdBiodata(0);
                cd(correctionData.info.correctFiles{i,1});
                
                %read correctionMovie
                corrMovieName = correctionData.info.correctFiles{i,2};
                if strcmp(corrMovieName(end),'d');
                    corrMov = r3dread(corrMovieName,[],[],MOVIERANGE,0);
                else
                    corrMov = readmat(corrMovieName);
                end
                %average over timepoints
                corrImgTmp = mean(corrMov,5);
                %delete correctionMove
                clear corrMov;
                %store corrImg (max 2 movies, so no problem with mean)
                corrImg = squeeze(mean(cat(5,corrImg,corrImgTmp),5));
                
            end
            %store only one slice
            backgroundCorrectionSlice = mean(corrImg,3);
            
            %return to movieDir
            cd(oldDir);
            
            % make backgroundimage
            backgroundCorrectionImage = repmat(backgroundCorrectionSlice,[1,1,numZ,1]);
    end
    
    %calculate minStart, maxEnd
    if ~isempty(correctionData)
        if ~isempty(correctionData.info.correctFrames)
            %some frames are not interesting
            minStart = correctionData.info.correctFrames(1) + 1;
            maxEnd   = numTimes - correctionData.info.correctFrames(2);
        else
            minStart = 1;
            maxEnd   = numTimes;
        end
    end
    
    
    %save correctionData if necessary
    if ~isempty(correctionData) & isfield(correctionData,'image')
        % do not save, because we have already calculated the image before!
    else
        
        
        %read header
        
        if isfield(correctionData,'header')
            r3dMovieHeader = correctionData.header;
            correctionData = rmfield(correctionData,'header');
        else
            r3dMovieHeader = readr3dheader(filename);
        end
        
        %save without header
        
        correctionData.image = backgroundCorrectionSlice;
        save('correctionData','correctionData');
        
        %correct time, movieLength if necessary
        if ~isempty(correctionData.info.correctFrames)
            % correct header.time. This is very likely not compatible with
            % multiple wavelengths.
            % Since it is a 1 by numzSlices*numTimepoints
            % vector, we have to calculate from where to where to take it
            firstNum = (correctionData.info.correctFrames(1))*numZ+1; 
            lastNum = (numTimes - sum(correctionData.info.correctFrames)) * numZ; %ms5+cF1+cF2-cF2
            r3dMovieHeader.Time = r3dMovieHeader.Time(firstNum:lastNum);
            r3dMovieHeader.numTimepoints = numTimes - sum(correctionData.info.correctFrames);
        end
        r3dMovieHeader.cropInfo = []; %for compatibility
        r3dMovieHeader.correctInfo = correctionData.info;
                
        %save movieHeader
        save('r3dMovieHeader','r3dMovieHeader');
    end
    
    
    
