function batchTest(selectROI,nFrames)

topDir=uigetdir(pwd,'Please select top directory containing all the project directories');
dirNames = dir(topDir);
dirNames(1:2)=[]; % get rid of . and .. directories in list

% get rid out of the list any element in dirNames that isn't a real
% directory (e.g. .txt file, .mat file)
countNonDirectories=zeros(length(dirNames)-2,1);
for i=1:length(dirNames)
    countNonDirectories(i)=isdir([topDir filesep dirNames(i).name]);
end

dirNames(countNonDirectories==0)=[];
nProj=length(dirNames);
for i=1:nProj
    runInfo.imDir=[topDir filesep dirNames(i).name filesep 'images'];
    runInfo.anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

    % kill old spot directory if it exists and make a new one
    spotDir=[runInfo.anDir filesep 'spot'];
    if isdir(spotDir)
        rmdir(spotDir,'s');
    end
    mkdir(spotDir);

    % let user select a ROI for each project
    if selectROI==1
        [listOfImages]=searchFiles('.tif',[],runInfo.imDir,0);
        fileNameIm=[char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
        img=double(imread(fileNameIm));
        runInfo.roiMask=roipoly(img./max(img(:)));
    end

    save([spotDir filesep 'runInfo'],'runInfo');

end

% run eb-detector and frame integrator
for i=1 %:nProj
    spotDir=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'spot'];
    temp=load([spotDir filesep 'runInfo.mat']);
    runInfo=temp.runInfo;
    
    if nargin<2
        [listOfImages]=searchFiles('.tif',[],runInfo.imDir,0);
        nFrames=size(listOfImages,1); % use all the frames
    elseif length(nFrames)==nProj
        nFrames=nFrames(i); % use value i from nFrames-vector
    elseif length(nFrames)==1
        % nothing, do same number of frames for all projects
    else
        error('nFrames should be one number or nProjects-vector or left out');
    end
    
    
    % right now debug is basically meaningless...
    eb1SpotDetector(runInfo,nFrames,0);
    
    % run frame integrator
    imIntegrator(runInfo,nFrames);
    
    close all
    
    % get first integrated image (over all frames) for overlay
    [listOfIntImages]=searchFiles('meanImg',[],[spotDir filesep 'intIm'],0);
    fileNameIm=[char(listOfIntImages(1,2)) filesep char(listOfIntImages(1,1))];
    img=double(imread(fileNameIm));
    
    
    saveResults.dir=spotDir;
    temp=load([spotDir filesep 'movieInfo.mat']);
    movieInfo=temp.movieInfo;
    % run tracking
    scriptTrackGeneral
    
    % plot tracks in magenta, +:trackStart, square:trackEnd
    plotTracks2D(tracksFinal,[],'m',[],1,1,img,0,0);

    % plot detected centroids on top using jet colormap to distinguish
    % frames
    spotsOverTime(runInfo,nFrames);
    
    saveas(gcf,[spotDir filesep 'forwardTrackResults'])
    
    
    % plot timecourse of spots over integrated frames
    spotsOverTime(runInfo,nFrames)

    close all
end



%[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
% nImages=25; img=0;
% for i=1:nImages
%     fileNameIm=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
%     if i==1
%         singleImg=Gauss2DBorder(double(imread(fileNameIm)),1);
%     end
%     img=img+double(imread(fileNameIm));
% end
% 
% meanImg=img./nImages;
% imshow(singleImg-meanImg,[])

% img.data=meanImg;
% img.perm='C';
% [resp,ori,nse,scaleMap,maxMap]=imLineDetect(img,[3 6],1,.99);
