function [CCPstruct,nascentAdhInfo,focalAdhInfo] = getDistanceCCPToAdhesion(pathForTheMovieDataFile,bandwidth,iCCP,iPax)

if nargin==1
    bandwidth = 5;
end
if nargin<3
    iCCP = 1; %assumed
end
if nargin<3
    iPax=3;
end
% Load the MovieData
if isa(pathForTheMovieDataFile, 'MovieData')
    movieData = pathForTheMovieDataFile;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
    movieData = MovieData.load(movieDataPath,false);
end
% run detectMovieNascentAdhesion
plotGraph = false;
[nascentAdhInfo,focalAdhInfo] = detectMovieNascentAdhesion(movieData,bandwidth,iPax,plotGraph);

%% --------------- Set up ---------------%%
nChan=numel(movieData.channels_);
% psfSigma = 1.2; %hard coded for nascent adhesion
psfSigma = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, movieData.getChannel(iCCP).emissionWavelength_*1e-9);

% Set up the output directories
outputFilePath = [movieData.outputDirectory_ filesep 'CCP_distanceAnalysis'];
imgPath = [outputFilePath filesep 'imgs'];
dataPath = [outputFilePath filesep 'data'];
tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];

if ~exist(figPath,'dir')
    mkdir(imgPath);
    mkdir(dataPath);
    mkdir(tifPath);
    mkdir(figPath);
end

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting CCPs...')

roiMask = movieData.getROIMask;
CCPstruct(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[], 'numberCCP',[],'cellArea',[],'CCPdensity',[],'distToAdh',[],'closeToNA',[],'idNA',[],'idFA',[]);

h1=figure;
jformat = ['%.' '3' 'd'];
% Changed it for isometric detection for nascent adhesion detection
pixSize = movieData.pixelSize_;
% iCCP = 1;
disp(['CCP channel was assumed to be in channel ' num2str(iCCP) '.'])
disp('Results will be saved under:')
disp(outputFilePath);
    
for j=1:movieData.nFrames_
    Iccp=double(movieData.channels_(iCCP).loadImage(j));
    maskProc = movieData.getProcess(movieData.getProcessIndex('MaskRefinementProcess'));
    % if there are masks for more than one channels, combine them.
    if length(maskProc.checkChannelOutput)>1
        %Combine the the multiple masks to one
        maskEach = arrayfun(@(x) maskProc.loadChannelOutput(x,j),find(maskProc.checkChannelOutput),'UniformOutput',false);
        maskAll=reshape(cell2mat(maskEach),size(Iccp,1),size(Iccp,2),[]);
        mask = any(maskAll,3);
    elseif length(iChans)==1
        mask = maskProc.loadChannelOutput(iPax,j); % 1 is CCP channel
    end
    ultimateMask = roiMask(:,:,j) & mask;
    pstruct = pointSourceDetection(Iccp, psfSigma, 'Alpha',0.05,'Mask', ultimateMask);
    % filter out points where overall intensity is in the noise level
    pixelIntenMargin = Iccp(~ultimateMask);
    maxIntBg=quantile(pixelIntenMargin,0.99);
    psInt = pstruct.A+pstruct.c;
    idxSigCCP = psInt>maxIntBg;
    
    CCPstruct(j).xCoord = [round(pstruct.x(idxSigCCP)'), round(pstruct.x_pstd(idxSigCCP)')];
    CCPstruct(j).yCoord = [round(pstruct.y(idxSigCCP)'), round(pstruct.y_pstd(idxSigCCP)')];
    CCPstruct(j).amp = [pstruct.A(idxSigCCP)', pstruct.A_pstd(idxSigCCP)'];
    CCPstruct(j).numberCCP = length(pstruct.x(idxSigCCP));
    CCPstruct(j).cellArea = sum(ultimateMask(:))*(pixSize/1000)^2; % in um^2
    CCPstruct(j).CCPdensity = CCPstruct(j).numberCCP/CCPstruct(j).cellArea; % number per um2

    % Closest distance analysis: distance to NAs
    [idxNAfromCCP, distToNA] = KDTreeClosestPoint([nascentAdhInfo(j).xCoord(:,1), nascentAdhInfo(j).yCoord(:,1)],...
                                                    [CCPstruct(j).xCoord(:,1),CCPstruct(j).yCoord(:,1)]);
    CCPstruct(j).distToAdh = distToNA;
    CCPstruct(j).closeToNA = true(size(distToNA)); %true means it is close to NA than FA (false: close to FA more than NA)
    CCPstruct(j).idNA = idxNAfromCCP;
    
    % distance to FAs
    allBoundFAs=cell2mat(focalAdhInfo(j).boundFA);
    integratedFApoints = [focalAdhInfo(j).xCoord(:,1), focalAdhInfo(j).yCoord(:,1);...
                                                                             allBoundFAs(:,2), allBoundFAs(:,1)];
    [idxFAfromCCP, distToFA] = KDTreeClosestPoint(integratedFApoints,...
                                                    [CCPstruct(j).xCoord(:,1),CCPstruct(j).yCoord(:,1)]);
    CCPstruct(j).distToAdh(distToFA<distToNA) = distToFA(distToFA<distToNA);
    CCPstruct(j).closeToNA(distToFA<distToNA) = false; %(false: close to FA more than NA)
    CCPstruct(j).idFA = idxFAfromCCP;
    
    % plotting detected adhesions
    Ipax=double(movieData.channels_(iPax).loadImage(j));
    dIpax =double(Ipax)/max(max(Ipax)); 
    dIccp = double(Iccp)/max(max(Iccp));
    combI(:,:,1) = dIpax;
    combI(:,:,2) = dIccp;
    combI(:,:,3) = zeros(size(dIpax));
    imshow(combI,[]); hold on
    plot(CCPstruct(j).xCoord(:,1),CCPstruct(j).yCoord(:,1),'go')
    plot(nascentAdhInfo.xCoord(:,1),nascentAdhInfo.yCoord(:,1),'y.')
    nFA = focalAdhInfo.numberFA;
    for k = 1:nFA
        adhBoundary = focalAdhInfo.boundFA{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'Color','r', 'LineWidth', 0.5) %adhesion boundary
    end
    
    plot([CCPstruct(j).xCoord(CCPstruct(j).closeToNA,1) nascentAdhInfo(j).xCoord(CCPstruct(j).idNA(CCPstruct(j).closeToNA),1)]',...
        [CCPstruct(j).yCoord(CCPstruct(j).closeToNA,1) nascentAdhInfo(j).yCoord(CCPstruct(j).idNA(CCPstruct(j).closeToNA),1)]','y')
    plot([CCPstruct(j).xCoord(~CCPstruct(j).closeToNA,1) integratedFApoints(CCPstruct(j).idFA(~CCPstruct(j).closeToNA),1)]',...
        [CCPstruct(j).yCoord(~CCPstruct(j).closeToNA,1) integratedFApoints(CCPstruct(j).idFA(~CCPstruct(j).closeToNA),2)]','r')
    hgexport(h1,strcat(tifPath,'/imgCCPtoAdh',num2str(j,jformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h1,strcat(figPath,'/imgCCPtoAdh',num2str(j,jformat)),'-v7.3')
    hold off
    close(h1)
end
save([dataPath filesep 'AdhInfo.mat'],'nascentAdhInfo','focalAdhInfo','CCPstruct','-v7.3');

disp('Finished detecting and linking objects...')