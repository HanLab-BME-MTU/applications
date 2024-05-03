function MO = createBestFocusedImageBF(pathToNDFile,iChan)
%function MO = createBestFocusedImageBF(pathToNDFile,iChan) accepts
%pathToNDFile, import it using bfImport, find the best focused-image using
%function findBestFocusFromStack, and produce MovieObject MO that has
%channels with only best-focused images.
% input
%   pathToNDFile        a character path toward a nd file
%   iChan               the channel number for Bead channel
% output
%   MO                  movieData or movieList  with channels with the
%                       focused images
%   Sangyoon Han, Nov 2021
% e.g.: MO = createBestFocusedImageBF
%% Read the nd file
% if it is not entered as an input argument, ask the user.
if nargin<1
    [fileSFolders, pathSFolder] = uigetfile('*.nd','Select a nd file.');
    pathToNDFile = [pathSFolder filesep fileSFolders];
else
    [pathSFolder,name,ext] = fileparts(pathToNDFile);
%     fileSFolders=[name ext];
end
MD = bfImport(pathToNDFile);

if nargin<2
    nChan = numel(MD(1).channels_);
    if nChan >1
        % Show each channgel
        hf = figure;
        for ii=1:nChan
            subplot(1,nChan,ii)
            imshow(MD(1).channels_(ii).loadImage(1),[])
            title(['Channel ' num2str(ii)])
        end
        % Ask user for channel index
        iChan = input('Enter bead channel number (default: the final one): ');
        if isempty(iChan)
            iChan = numel(refMD.channels_);
        end
        close(hf)
    elseif nChan==1
        iChan=1;
    end
end
% Going over each movie
nMovies = numel(MD);
applySobel = true;
for jj=1:nMovies
    curMD = MD(jj);
    % 
    if curMD.zSize_>1
        curBeadChan = curMD.channels_(iChan);
        % store it somewhere
        curPath = curMD.outputDirectory_;
        % Make the channel folder
        curNChan = numel(curMD.channels_);
        curChanAll = cell(curNChan,1);
        for kk=1:curNChan
            curChanPath = [curPath filesep 'Chan' num2str(kk)];
            curChanAll{kk} = curChanPath;
            mkdir(curChanPath)
        end
        thresVariance=0.8;
        numFrames = curMD.nFrames_;
        for ii=1:numFrames
            % find the best focus
            curBeadStack = curBeadChan.loadStack(ii);
            averagingRange = findBestFocusFromStack(curBeadStack,thresVariance,applySobel);
            %curBestImage = curBeadStack(:,:,median(averagingRange));
            %figure, imshow(curBestImage,[])
            for kk=1:curNChan
                if kk==iChan
                    curBeadStackChosen = curBeadStack(:,:,averagingRange);
                    curBestImage = mean(curBeadStackChosen,3); 
                else
                    curChan = curMD.channels_(kk);
                    curStack = curChan.loadStack(ii);
                    curStackChosen = curStack(:,:,averagingRange);
                    curBestImage = mean(curStackChosen,3); 
                end
                meanImgPath = [curChanAll{kk} filesep 'img' num2str(ii) '.tif'];
                imwrite(uint16(curBestImage),meanImgPath,'Compression','none')
            end
        end
        % Registering it to a new MD file
        %channel(curNChan) = Channel();
        for kk=1:curNChan
            curChan = curMD.channels_(kk);
            channel(kk) = Channel(curChanAll{kk});
            try
                channel(kk).fluorophore_=curChan.name_;
                channel(kk).emissionWavelength_=name2wavelength(curChan.name_)*1e9;
                channel(kk).imageType_='Confocal'; % this should be changed later.
            catch
                disp(['Channel ' num2str(kk) ' is ' curChan.name_])
            end
        end
        MDnew(jj) = MovieData(channel,curPath);

        % Set the path where to store the MovieData object.
        MDnew(jj).setPath(curPath);
        MDnew(jj).setFilename('movieData.mat');

        % Set some additional movie properties
        MDnew(jj).numAperture_= curMD.numAperture_; %1.49;
        if isempty(curMD.pixelSize_)
            magObj = input('What is the objective magnification (e.g. 20)?: ');
            MDnew(jj).pixelSize_= 6500/magObj;%71;
        else
            MDnew(jj).pixelSize_= curMD.pixelSize_;%71;
        end
        MDnew(jj).camBitdepth_=curMD.camBitdepth_;
        if isempty(curMD.timeInterval_)
            tInt = input('What is the time interval in seconds?: ');
            MDnew(jj).timeInterval_= tInt;%71;
        else
            MDnew(jj).timeInterval_= curMD.timeInterval_;
        end
        MDnew(jj).notes_= 'Movie with only focused images '; 

        % Run sanityCheck on MovieData. 
        % Check image size and number of frames are consistent. 
        % Save the movie if successfull
        MDnew(jj).sanityCheck;
    else
        disp('This MD was not z-stack. You can use this directly after checking focused images')
    end
end
if nMovies==1
    MO= MDnew;
else
    MO = MovieList(MDnew,pathSFolder);
    MO.movieListFileName_ = 'movieList.mat';
    MO.sanityCheck;
end
%     if numel(averagingRange)>5
%         applySobel=false;
%         averagingRange = findBestFocusFromStack(curBeadStack,thresVariance,applySobel);
%         if numel(averagingRange)>5
%             averagingRange = findBestFocusFromStack(curBeadStack,thresVariance,applySobel,'amp');
%         end
%     end
% else
%     meanRefImg = curBeadChan.loadImage(1);
% end
end

