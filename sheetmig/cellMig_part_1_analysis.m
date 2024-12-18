function []=cellMig_part_1_analysis(r,rMin,rMax,smoothFac,imageFileListNuclei,imageFileListPhase,showMovie,sglCell,doNuclei)

if nargin<1
    r   =input('Specify the minimal nuclei radius    r=[ 5pix]: ');
    rMin=input('Specify the minimal search radius rMin=[    r]: ');
    rMax=input('Specify the maximal search radius rMax=[50pix]: ');
    smoothFac=input('Specify smoothing factor for smoother edges [2]: ');

    if isempty(r)
        r=5;
    end
    if isempty(rMin)
        rMin=r;
    end
    if isempty(rMax)
        rMax=50;
    end
    if isempty(smoothFac)
        smoothFac=2;
    end
end

resultDir=pwd;
if ~isdir(resultDir)
    mkdir(resultDir);
end

if nargin<9 || isempty(doNuclei)
    doNuclei=1;
end

if doNuclei
    if nargin < 5 || isempty(imageFileListNuclei)
        try
            imageFileListNuclei=[pwd,filesep,'Nuclei'];
            imageFileListNuclei=getFileListFromFolder(imageFileListNuclei);
        catch exception
            [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
                'Select First nuclei Image');
            
            if ~ischar(filename) || ~ischar(pathname)
                return;
            end
            
            imageFileListNuclei = getFileStackNames([pathname filesep filename]);
        end
        
    elseif isdir(imageFileListNuclei)
        imageFileListNuclei=getFileListFromFolder(imageFileListNuclei);
    else
        isValid = 1;
        for frame = 1:numel(imageFileListNuclei)
            isValid = isValid && exist(imageFileListNuclei{frame}, 'file');
        end
        if ~isValid
            error('Invalid input files.');
        end
    end
    
    [~, ~, fno , ~]=getFilenameBody(imageFileListNuclei{1});
    if str2double(fno)~=1
        display('Collapsed file stack! First frame had frame no > 1!')
        imageFileListNuclei=collapseFileStack(imageFileListNuclei, -1);
    end
end


if nargin<6 || isempty(imageFileListPhase)
    imageFileListPhase=[pwd,filesep,'Phase'];
    try
        imageFileListPhase=getFileListFromFolder(imageFileListPhase);
    catch exception2
        [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
            'Select First phase Image');
        
        if ~ischar(filename) || ~ischar(pathname)
            return;
        end
        
        imageFileListPhase = getFileStackNames([pathname filesep filename]);
    end
elseif isdir(imageFileListPhase)
    imageFileListPhase=getFileListFromFolder(imageFileListPhase);
else
    isValid = 1;
    for frame = 1:numel(imageFileListPhase)
        isValid = isValid && exist(imageFileListPhase{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end
[~, ~, fno , ~]=getFilenameBody(imageFileListPhase{1});
if str2double(fno)~=1
    display('Collapsed file stack! First frame had frame no > 1!')
    imageFileListPhase=collapseFileStack(imageFileListPhase, -1);
end

if nargin<8 || isempty(sglCell)
    sglCell=0;
end


if doNuclei
    % find nuclei:
    display('Detect nuclei edge:...')
    [nuclei,movieInfo,dPix]=nucleiDetect(imageFileListNuclei,r,1,'sobel',[],0.01);
    display('Done! Save results:...')
    save('xDetectNuclei.mat','nuclei','movieInfo','dPix','imageFileListNuclei','imageFileListPhase');
    close all;
else
    dPix=[];
end


% find the wound edge
%!!! To be done: There might be holes at the boundary that have to be
%filled!
display('Detect wound edge:...')
if sglCell
    [sheetMask,sheetBnD,sheetEdge,~,~,toDoList]=createSglCellMask(imageFileListPhase,r,dPix);    
else
    [sheetMask,sheetBnD,sheetEdge,~,~,toDoList]=woundEdgeDetect(imageFileListPhase,r,smoothFac,'canny',2,dPix,'noholes');
end
display('Done! Save results:...')
save('xDetectEdge.mat','sheetMask','sheetBnD','sheetEdge','toDoList','sglCell','smoothFac','-v7.3');
% The fields: ,'cellDistFromEdge','distGrad' are not saved!
close all;

if doNuclei
    % before tracking crop the nuclei that are within the cell sheet:
    display('Cut off out of sheet nuclei:...')
    [nucleiCrop,movieInfoCrop]=cropNucleiWithMask(imageFileListNuclei,r,nuclei,sheetMask,toDoList);
    display('Done! Save results:...')
    save('xDetectNuclei.mat','nuclei','movieInfo','dPix','imageFileListNuclei','imageFileListPhase','nucleiCrop','movieInfoCrop','-v7.3');
    
    % track nuclei:
    display('Track nuclei:...')
    %[tracksFinal]=scriptTrackNuclei(movieInfoCrop,rMin,rMax,resultDir);
    [tracksFinal]=scriptTrackNucleiWithGapCl(movieInfoCrop,rMin,rMax,resultDir);
    display('Done! Save results:...')
    close all;
    
    if showMovie==1
        % show and save movie:
        movieFileName='mov_trackNuclei.mov';
        display('Create tracker movie:...')
        overlayTracksMovieNew(tracksFinal,[],10,1,movieFileName,[],1,0,0,[],0,1,[],0,0,imageFileListNuclei)
        % move the movie into the current dir:
        if isdir('Nuclei')
            try
                movefile(['Nuclei',filesep,movieFileName],movieFileName);
            catch exception
                display('Couldnt find movie file, you have to find it yourself');
            end
        end
    end
end
save('xParameters.mat','r','rMin','rMax','smoothFac');
display('All done!')
