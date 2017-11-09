function [maskList,mask,fixedFrames] = calcCellBoundary(imageDirectory,outputDirectory,doPlot,isMultiChannel,fixFrameUp, fixFrameDown)
%CALCCELLBOUNDARY applies mutli-otsu threshold to determine cell boundary
% Mask is best suited for images in which cell is large relative to image
% dimensions
% Synopsis: [maskList,mask] = calcCellBoundary(image)
% Input:
%   firstImageFile - Name, including full path, of first image tiff file.
%    Ex: '/home2/avega/Desktop/CD36 Fyn/+TSP-1/image_0001.tif'
%    Optional. User will be prompted to choose file if not supplied.
%
%   outputDirectory - Path of directory where results will be saved
%
%   doPlot - 1 to show result, 0 otherwise. Default: 0.
%
%   isMultiChannel - If any value is given, will prompt user to input
%   channel where image is found
%
%   fixFrameUp - Input array containing the number of each frame which
%   requires higher threshold
%
%   fixFrameDown - Input array containing the number of each frame which
%   requires lower threshold
%
% Output:
%   maskList- cell array in which each cell contains n x 2 list of all 
%   pixels within cellBoundary for a given frame
%
%   mask - three dimensional matrix where each slice in z-dimension
%   contains mask for corresponding frame
%
%   fixedFrames- structure containing numbers of frame which required
%   lower or higher threshold for later reference
%% Input

if nargin < 1 || isempty(imageDirectory)
    [fName,dirName] = uigetfile('*.tif','PLEASE SPECIFY FIRST IMAGE IN STACK');
else
    fName = imDir(imageDirectory);
%     if iscell(firstImageFile)
%         [fpath,fname,fno,fext]=getFilenameBody(firstImageFile{1});
%         dirName=[fpath,filesep];
%         fName=[fname,fno,fext];
%     elseif ischar(firstImageFile)
%         [fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
%         dirName=[fpath,filesep];
%         fName=[fname,fno,fext];
%     end
end
outFileListCnt = vertcat(fName.name);
numFiles = length(outFileListCnt);

if nargin < 3 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 4 || isempty(isMultiChannel)
    channel = 1;
else
    
    channel = isMultiChannel;
end

% Initialize Movie
 writerObj = VideoWriter(outputDirectory);
 writerObj.FrameRate=1;
 open(writerObj);

 maskList = cell(numFiles,1);
 fixedFrames = struct('HigherThreshold',fixFrameUp,'LowerThreshold',fixFrameDown);
%% Masking Process
for a = 1:numFiles
    image  = imread(outFileListCnt(a,:),channel);
    %keep original image for plotting
    image0 = image;

    % Apply wide gaussian filter
     image = filterGauss2D(image,3);
%     [image, ~, ~, ~] = steerableDetector(double(image), 2, 2);
     
     
     
     % Compute the thresholds
    Nvals = [1:20];
    metric = zeros(length(Nvals),1);
    for i = 1:length(Nvals)
%         if a == 14 && i == 2
%             metric(i) = metric(i-1);
%             continue
%         end
        [~, metric(i)] = multithresh(image, Nvals(i) );
%         disp(['N = ' int2str(Nvals(i)) '  |  metric = ' num2str(metric(i))]);
    end

     
    %Apply multi-Otsu threshold on image
    thresh = multithresh(image,Nvals(find(metric == (max(metric)))));
    %Attempt to find largest gap in thresholds
    diffThresh = zeros(length(thresh)-1,1);
    for i = 1:(length(thresh)-1)
       diffThresh(i) = thresh(i+1)/thresh(i);
    end
    threshLevel =1 + (find(diffThresh == max(diffThresh)));
    
     if isempty(find(fixFrameUp == a)) && isempty(find(fixFrameDown == a))
         test = image >= thresh(threshLevel);
     elseif isempty(find(fixFrameDown == a))
         test = image >= thresh(threshLevel+1);
     else
         test = image >= thresh(threshLevel-1);
     end
    
    test = imfill(test,'holes');
    

    % Find Connected Components and their areas
    [L, ~] = bwlabel(test, 8);
    STATS = regionprops(L, 'Area');


    %Find all areas and then only get indices of components that have at least
    %some fraction of the largest area (default: 0.5)
    areaInfo = cat(1,STATS.Area);
    sThreshold = find(areaInfo >= (max(areaInfo)/2));


    % Create new matrix with only the selected components, there's probably an
    % easier way of doing this...
    for i =1:length(sThreshold)
        maskTemp(:,:,i) = (L ==sThreshold(i));
        maskTemp(:,:,1) = maskTemp(:,:,1)+maskTemp(:,:,i);
    end

    [i,j] = find(maskTemp(:,:,1));
    maskListTemp(:,1) = i;
    maskListTemp(:,2) = j;
    
    maskList{a} = maskListTemp; 
    mask(:,:,a) = maskTemp(:,:,1); 
    if doPlot
        h = figure;
        imshow(image0,[]);
        hold on;
        maskBounds = bwboundaries(maskTemp(:,:,1));
        cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
    end
    % Make frame of movie from output
    scrsz = get(0,'ScreenSize');
    set(h,'Position',scrsz);
    testMovie  = getframe(h);
    writeVideo(writerObj,testMovie);
    close(h);
    clear maskListTemp maskTemp
end
close(writerObj);

save('/project/biophysics/jaqaman_lab/vegf_tsp1/touretLab/CtxB-CD36-Actin/maskingFileNoTSP.mat','mask','maskList','fixedFrames');
end