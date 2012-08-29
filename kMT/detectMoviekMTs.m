function detectMoviekMTs(movieData,varargin)

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

iProc=movieData.getProcessIndex('KMTDetectionProcess',1,0);
if isempty(iProc)
    movieData.addProcess(KMTDetectionProcess(movieData));
    iProc=numel(movieData.processes_);
end
detProc=movieData.processes_{iProc};


%Parse input, store in parameter structure
p = parseProcessParams(detProc,paramsIn);

% gfpIndx = strcmpi('gfp',{movieData.channels_.fluorophore_});

% Detectino mask parameters
sigmaY=1.5;
sigmaX=2*sigmaY;
% sigmaY = movieData.channels_(gfpIndx).psfSigma_;
% sigmaX=2*sigmaY;
kSigma=2;


% Get groupign proces index
iGroupProc = movieData.getProcessIndex('SisterGroupingProcess',Inf,0);
assert(~isempty(iGroupProc));
groupProc=movieData.processes_{iGroupProc};

groupChannelOutput = groupProc.checkChannelOutput;
assert(any(groupChannelOutput));
sisterList = groupProc.loadChannelOutput(find(groupChannelOutput,1));
nPairs=numel(sisterList);


% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = groupProc.outFilePaths_{1,find(groupChannelOutput,1)};
end
detProc.setInFilePaths(inFilePaths);

%Set up the mask directories as sub-directories of the output directory
outFilePaths= cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;
    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    outFilePaths{2,i} = [p.OutputDirectory filesep 'filtered_images_for_channel_' num2str(i)];

    %Check/create directory
    mkClrDir(outFilePaths{1,i});
    mkClrDir(outFilePaths{2,i})

end
detProc.setOutFilePaths(outFilePaths);

disp('Starting decting kMT mask');

% get difference of Gaussians image for each frame and standard deviation
% of the cell background, stored in stdList

nFrames=movieData.nFrames_;
maxIntensity=(2^movieData.camBitdepth_)-1;
roiMask=movieData.getROIMask();
stdList=NaN(nFrames,1);
logMsg=['Filtering images for comet detection for channel ' num2str(p.ChannelIndex)];
progressText(0,logMsg);

for iFrame=1:nFrames
    
    % Loading image
    im = movieData.channels_(p.ChannelIndex).loadImage(iFrame)/maxIntensity;
    
    % Combine region of interest with preselected segmentation output
    mask=roiMask(:,:,iFrame);
    
    
    % Difference of Gaussians
    filterDiff=filterGauss2D(im,p.sigma1)-filterGauss2D(im,p.sigma2);
    filterDiff(~mask)=NaN;
    
    % Save filtered images on disk (avoid memory errors)
    save(fullfile(outFilePaths{2,p.ChannelIndex},['filterDiff_' num2str(iFrame) '.mat']),'filterDiff');
    
    % Perform maximum filter and mask out significant pixels
    thFilterDiff = ordfilt2(filterDiff,9,ones(3,3));
    threshold = thresholdOtsu(thFilterDiff)/3 + ...
        thresholdRosin(thFilterDiff)*2/3;
    stdList(iFrame)=nanstd(filterDiff(thFilterDiff<threshold));
    %         stdList(i)=nanstd(filterDiff(:));
    
    % Update progress status
    progressText(iFrame/nFrames,'Filtering images for comet detection');
%     if ishandle(wtBar) && mod(i,5)==0, waitbar(frac,wtBar,logMsg);  end
end

meanStd = arrayfun(@(x) mean(stdList(max(1,x-1):min(nFrames,x+1))),1:nFrames);
logMsg='Detecting comets';
progressText(0,logMsg);

for iPair=1:numel(sisterList)
    sisterList(iPair).kMTcoords1=NaN(nFrames,3);
    sisterList(iPair).kMTamp1=NaN(nFrames,3);
    sisterList(iPair).kMTcoords2=NaN(nFrames,3);
    sisterList(iPair).kMTamp2=NaN(nFrames,3);

end

for iFrame=1:nFrames
    % loop thru frames and detect
    
    % Reload band-pass filtered images
    s=load(fullfile(outFilePaths{2,p.ChannelIndex},['filterDiff_' num2str(iFrame) '.mat']));
    filterDiff=s.filterDiff;
    stepSize=p.multFactorStepSize*meanStd(iFrame);
    thresh= p.multFactorThresh*meanStd(iFrame);
    
    mask=false(movieData.imSize_);
    validPairs =  ~arrayfun(@(x) all(isnan(x.coords1(iFrame,:))),sisterList);

    for iPair=find(validPairs)'
        
        
        dP = sisterList(iPair).sisterVectors(iFrame,1:2);
        pos1 = sisterList(iPair).coords1(iFrame,1:2);
        pos2 = sisterList(iPair).coords2(iFrame,1:2);
        e = dP/norm(dP);
        theta = atan2(dP(2),dP(1));
        [xRange,yRange,nzIdx] = anisoGaussian2DSupport(pos1(1)-e(1)*kSigma*sigmaX/2,...
            pos1(2)-e(2)*kSigma*sigmaX/2,sigmaX,sigmaY,theta,kSigma,movieData.imSize_);
        %         plot([min(xRange) max(xRange) max(xRange) min(xRange) min(xRange)],...
        %             [min(yRange) min(yRange) max(yRange) max(yRange) min(yRange)],'-k');
        [X,Y]=meshgrid(xRange,yRange);
        ind=sub2ind(movieData.imSize_,Y(nzIdx),X(nzIdx));
        mask(ind)=true;
        movieInfo =detectComets(filterDiff.*mask,stepSize,thresh);
        if isempty(movieInfo.xCoord)
            sisterList(iPair).kMTcoords1(iFrame,1:2)=[NaN NaN];
        else
            kMTcoords = [movieInfo.xCoord(:,1) movieInfo.yCoord(:,1)];
            dl = kMTcoords(:,1:2) - repmat(pos1,[size(kMTcoords,1) 1]);
            distances = sqrt(sum(dl.^2,2));
            [~,index]=min(distances);
            sisterList(iPair).kMTcoords1(iFrame,:)=[kMTcoords(index,:) 0];
            sisterList(iPair).kMTamp1(iFrame,1)=movieInfo.int(index);

        end
            
        [xRange,yRange,nzIdx] = anisoGaussian2DSupport(pos2(1)+e(1)*kSigma*sigmaX/2,...
            pos2(2)+kSigma*e(2)*sigmaX/4,sigmaX,sigmaY,theta,kSigma,movieData.imSize_);
        %         plot([min(xRange) max(xRange) max(xRange) min(xRange) min(xRange)],...
        %             [min(yRange) min(yRange) max(yRange) max(yRange) min(yRange)],'-k')
        
        [X,Y]=meshgrid(xRange,yRange);
        ind=sub2ind(movieData.imSize_,Y(nzIdx),X(nzIdx));
        mask(ind)=true;
        
        movieInfo =detectComets(filterDiff.*mask,stepSize,thresh);
        if isempty(movieInfo.xCoord)
            sisterList(iPair).kMTcoords2(iFrame,1:2)=[NaN NaN];
        else
            kMTcoords = [movieInfo.xCoord(:,1) movieInfo.yCoord(:,1)];
            dl = kMTcoords(:,1:2) - repmat(pos2,[size(kMTcoords,1) 1]);
            distances = sqrt(sum(dl.^2,2));
            [~,index]=min(distances);
            sisterList(iPair).kMTcoords2(iFrame,:)=[kMTcoords(index,:) 0];
            sisterList(iPair).kMTamp2(iFrame,1)=movieInfo.int(index);

        end
    end
    
  
%     % Detect comets using watershed detection
%     movieInfo(iFrame,1) = detectComets(filterDiff.*mask,stepSize,thresh);
%     
    % Update progress status
    progressText(iFrame/nFrames,logMsg);
    
%     imwrite(mask,fullfile(outFilePaths{1,p.ChannelIndex},['mask_' num2str(iFrame) '.tif']));

end

% movieInfo=allMovieInfo(:,p.ChannelIndex); %#ok<NASGU>
save(outFilePaths{1,p.ChannelIndex} ,'sisterList','stdList');
disp('Finished creating kMT mask');