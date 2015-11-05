function [] = makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmax,nWorkers)
%function [] = makeMoviesAdhesionTracks(tracks,pathForColocalization,tmax) makes
% tif series for movie with classified adhesion tracks
close all
disp('Loading raw files ...')
band = 10;
imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
imgMap = imgMap.paxImgStack;
tMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
tMap = tMap.tMap;
cropInfo =  load([pathForColocalization filesep 'data'  filesep 'cropInfo.mat'],'cropInfo');
cropInfo = cropInfo.cropInfo;
numFrames = size(imgMap,3);
width_b = (cropInfo(3)-cropInfo(1)-2*band+1); %width with bands cut
height_b = cropInfo(4)-cropInfo(2)-2*band+1; %height with bands cut
disp('Predicting with classifier ...')
tic
T =  load([pathForColocalization filesep 'data'  filesep 'trainedClassifier.mat'],'trainedClassifier');
trainedClassifier = T.trainedClassifier;
toc
% trainedClassifier = trainClassifierNA(T);
[~,allData] = extractFeatureNA(tracksNA);
allDataClass = predict(trainedClassifier,allData);

moviePath = [pathForColocalization filesep 'movies'];
movieEpsFluoPath = [moviePath filesep 'epsFluo'];
movieTifFluoPath = [moviePath filesep 'tifFluo'];
movieFigFluoPath = [moviePath filesep 'figFluo'];
movieEpsForcePath = [moviePath filesep 'epsForce'];
movieTifForcePath = [moviePath filesep 'tifForce'];
movieFigForcePath = [moviePath filesep 'figForce'];
movieRawFluoPath = [moviePath filesep 'rawFluo'];
movieRawForcePath = [moviePath filesep 'rawForce'];
if ~exist(moviePath,'dir') || ~exist(movieRawForcePath,'dir') 
    mkdir(moviePath);
    mkdir(movieEpsFluoPath);
    mkdir(movieTifFluoPath);
    mkdir(movieFigFluoPath);
    mkdir(movieEpsForcePath);
    mkdir(movieTifForcePath);
    mkdir(movieFigForcePath);
    mkdir(movieRawFluoPath);
    mkdir(movieRawForcePath);
end
iiformat = ['%.' '3' 'd'];
if nargin<3 || isempty(tmax)
    tmax = 0.8*max(tMap(:));
end
tmin = min(tMap(:));
% h1 = figure('color','w'); % for fluorescence
% set(h1, 'Position', [100 50 width_b height_b])
% axIMap = axes;
% h2=figure('color','w'); % for force
% set(h2, 'Position', [100+width_b 50 width_b*1.25 height_b])
%% testing spmd
% spmd
% %     N = 10;
% %     X = magic(N,N,N);          % Replicated on every worker
% %     C1 = codistributed(X); % Partitioned among the workers
% %     L= getLocalPart(C1)
%     results = zeros(1, numDataSets, codistributor()); 
%     for i = drange(1:numDataSets)
%         load(['\\central\myData\dataSet' int2str(i) '.mat'])
%         results(i) = processDataSet(i); 
%     end 
%     res = gather(results, 1); 
%     if labindex == 1
%         plot(1:numDataSets, res);
%         print -dtiff -r300 fig.tiff;
%         save \\central\myResults\today.mat res
%     end
% end
%% Saving imgMap and tMap in its single frames for parallel processing
for ii=1:numFrames
    curImgFrame = imgMap(:,:,ii);
    save(strcat(movieRawFluoPath,'/flscRaw',num2str(ii,iiformat),'.mat'),'curImgFrame')
    curForceFrame = tMap(:,:,ii);
    save(strcat(movieRawForcePath,'/forceRaw',num2str(ii,iiformat),'.mat'),'curForceFrame')
end
%% Looping
% for ii=1:numFrames
% parpool(100)
% imgMapDist = distributed(imgMap);
% tMapDist = distributed(tMap);
% numFrames = 3;
if nargin<4
    nWorkers=100;
end
if isempty(gcp('nocreate'))
    parpool(nWorkers)
end
spmd (nWorkers)
    for k=1:floor(numFrames/nWorkers)
        h1 = figure('color','w'); % for fluorescence
        axIMap = axes;
        hold off
        ii = (k-1)*nWorkers+labindex;
        curImgFrame=load(strcat(movieRawFluoPath,'/flscRaw',num2str(ii,iiformat),'.mat'));
        curImgFrame = curImgFrame.curImgFrame;
        imshow(curImgFrame,[],'Parent',axIMap); hold on
        set(h1, 'Position', [100 50 width_b height_b])
        set(axIMap,'XLim',[cropInfo(1)+band cropInfo(3)-band],'YLim',[cropInfo(2)+band cropInfo(4)-band])
        [htrackLine, htrackCircles] =drawClassifiedTracks(allDataClass,tracksNA,ii,axIMap);
        drawnow
        print('-depsc2', '-r150', strcat(movieEpsFluoPath,'/flscInten',num2str(ii,iiformat),'.eps'));
        print('-dtiff', '-r150', strcat(movieTifFluoPath,'/flscInten',num2str(ii,iiformat),'.tif'));
    %     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
%         hgsave(h2,strcat(movieFigFluoPath,'/flscInten',num2str(ii,iiformat)),'-v7.3')
        hold off

%         figure(h2), 
        h2=figure('color','w'); % for force
        set(h2, 'Position', [100+width_b 50 width_b*1.25 height_b])
        axTMap = subplot('Position',[0 0 0.8 1]); hold off
        curForceFrame=load(strcat(movieRawForcePath,'/forceRaw',num2str(ii,iiformat),'.mat'));
        curForceFrame = curForceFrame.curForceFrame;
        imshow(curForceFrame,[tmin tmax]); colormap jet; hold on;
        set(axTMap,'XLim',[cropInfo(1)+band cropInfo(3)-band],'YLim',[cropInfo(2)+band cropInfo(4)-band])
        for kk=1:numel(htrackLine)
            numCurLine = numel(htrackLine{kk});
            for kkk=1:numCurLine
                copyobj(htrackLine{kk}{kkk},axTMap)
            end
            %drawing
%             cellfun(@(x) copyobj(x,axTMap),htrackLine{kk})
        end
        for kk=1:numel(htrackCircles)
            numCurCircle = numel(htrackCircles{kk});
            for kkk=1:numCurCircle
                copyobj(htrackCircles{kk}{kkk},axTMap)
            end
%             cellfun(@(x) copyobj(x,axTMap),htrackCircles{kk})
        end
        subplot('Position',[0.8 0.1 0.1 0.8])
        axis tight
        caxis([tmin tmax]), axis off
        hc = colorbar('West');
        set(hc,'Fontsize',12)
        drawnow
        print('-depsc2', '-r150', strcat(movieEpsForcePath,'/force',num2str(ii,iiformat),'.eps'));
        print('-dtiff', '-r150', strcat(movieTifForcePath,'/force',num2str(ii,iiformat),'.tif'));
    %     hgexport(h1,strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
%         hgsave(h1,strcat(movieFigForcePath,'/force',num2str(ii,iiformat)),'-v7.3')
        hold off
    end
end
% close([h1 h2])
delete(gcp)
        