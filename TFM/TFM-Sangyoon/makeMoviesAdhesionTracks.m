function [] = makeMoviesAdhesionTracks(tracksNA,pathForColocalization,tmax,nWorkers,drawWithoutClasses)
%function [] = makeMoviesAdhesionTracks(tracks,pathForColocalization,tmax) makes
% tif series for movie with classified adhesion tracks
if nargin<5
%     drawWithoutClasses=false;
    drawWithoutClasses=10;
end
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
if drawWithoutClasses==10
    disp('Loading classifier ...')
    tic
    T =  load([pathForColocalization filesep 'data'  filesep 'trainedClassifier.mat'],'trainedClassifier');
    trainedClassifier = T.trainedClassifier;
    toc
% trainedClassifier = trainClassifierNA(T);
    disp('Predicting with classifier ...')
    [~,allData] = extractFeatureNA(tracksNA);
    allDataClass = predict(trainedClassifier,allData);
elseif drawWithoutClasses<10 %&& drawWithoutClasses>0
    allDataClass = cell(numel(tracksNA),1);
end
% if drawWithoutClasses
%     moviePath = [pathForColocalization filesep 'movies_noClass'];
% else
%     moviePath = [pathForColocalization filesep 'movies'];
% end
switch drawWithoutClasses
    case 0
        moviePath = [pathForColocalization filesep 'movies_noClass'];
    case 10
        moviePath = [pathForColocalization filesep 'movies'];
    case 1
        moviePath = [pathForColocalization filesep 'moviesForClass1'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1filtered');
            curIdsClassified=curIdsClassified.idGroup1filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup1');
            curIdsClassified=curIdsClassified.idGroup1;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group1';
            else
                allDataClass{p}='Group10';
            end
        end
    case 2
        moviePath = [pathForColocalization filesep 'moviesForClass2'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup2filtered');
            curIdsClassified=curIdsClassified.idGroup2filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup2');
            curIdsClassified=curIdsClassified.idGroup2;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group2';
            else
                allDataClass{p}='Group10';
            end
        end
    case 3
        moviePath = [pathForColocalization filesep 'moviesForClass3'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup3filtered');
            curIdsClassified=curIdsClassified.idGroup3filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup3');
            curIdsClassified=curIdsClassified.idGroup3;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group3';
            else
                allDataClass{p}='Group10';
            end
        end
    case 4
        moviePath = [pathForColocalization filesep 'moviesForClass4'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup4filtered');
            curIdsClassified=curIdsClassified.idGroup4filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup4');
            curIdsClassified=curIdsClassified.idGroup4;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group4';
            else
                allDataClass{p}='Group10';
            end
        end
    case 5
        moviePath = [pathForColocalization filesep 'moviesForClass5'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup5filtered');
            curIdsClassified=curIdsClassified.idGroup5filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup5');
            curIdsClassified=curIdsClassified.idGroup5;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group5';
            else
                allDataClass{p}='Group10';
            end
        end
    case 6
        moviePath = [pathForColocalization filesep 'moviesForClass6'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup6filtered');
            curIdsClassified=curIdsClassified.idGroup6filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup6');
            curIdsClassified=curIdsClassified.idGroup6;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group6';
            else
                allDataClass{p}='Group10';
            end
        end
    case 7
        moviePath = [pathForColocalization filesep 'moviesForClass7'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup7filtered');
            curIdsClassified=curIdsClassified.idGroup7filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup7');
            curIdsClassified=curIdsClassified.idGroup7;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group7';
            else
                allDataClass{p}='Group10';
            end
        end
    case 8
        moviePath = [pathForColocalization filesep 'moviesForClass8'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup8filtered');
            curIdsClassified=curIdsClassified.idGroup8filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup8');
            curIdsClassified=curIdsClassified.idGroup8;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group8';
            else
                allDataClass{p}='Group10';
            end
        end
    case 9
        moviePath = [pathForColocalization filesep 'moviesForClass9'];
        try
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup9filtered');
            curIdsClassified=curIdsClassified.idGroup9filtered;
        catch
            curIdsClassified=load([pathForColocalization filesep 'data' filesep 'idsClassified.mat'],'idGroup9');
            curIdsClassified=curIdsClassified.idGroup9;
        end            
        for p=1:numel(tracksNA)
            if ismember(p,find(curIdsClassified)')
                allDataClass{p}='Group9';
            else
                allDataClass{p}='Group10';
            end
        end
end
movieEpsFluoPath = [moviePath filesep 'epsFluo'];
movieTifFluoPath = [moviePath filesep 'tifFluo'];
movieFigFluoPath = [moviePath filesep 'figFluo'];
movieEpsForcePath = [moviePath filesep 'epsForce'];
movieTifForcePath = [moviePath filesep 'tifForce'];
movieFigForcePath = [moviePath filesep 'figForce'];
movieRawFluoPath = [moviePath filesep 'rawFluo'];
movieRawForcePath = [moviePath filesep 'rawForce'];
if ~exist(moviePath,'dir') || ~exist(movieRawForcePath,'dir') 
    mkClrDir(moviePath);
    mkClrDir(movieEpsFluoPath);
    mkClrDir(movieTifFluoPath);
    mkClrDir(movieFigFluoPath);
    mkClrDir(movieEpsForcePath);
    mkClrDir(movieTifForcePath);
    mkClrDir(movieFigForcePath);
    mkClrDir(movieRawFluoPath);
    mkClrDir(movieRawForcePath);
end
iiformat = ['%.' '3' 'd'];
if nargin<3 || isempty(tmax)
    tmax = 0.8*max(tMap(:));
end
tmin = min(tMap(:));
markerSize = 4;
% h1 = figure('color','w'); % for fluorescence
% set(h1, 'Position', [100 50 width_b height_b])
% axIMap = axes;
% h2=figure('color','w'); % for force
% set(h2, 'Position', [100+width_b 50 width_b*1.25 height_b])
%% testing spmd
% spmd
%     N = 10;
%     X = magic(N,N,N);          % Replicated on every worker
%     C1 = codistributed(X); % Partitioned among the workers
%     L= getLocalPart(C1)
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
%% Test parfor - this works
% cropInfo1=cropInfo(1);
% cropInfo2=cropInfo(2);
% cropInfo3=cropInfo(3);
% cropInfo4=cropInfo(4);
% parfor k=1:floor(numFrames)
%     figure('color','w'); % for fluorescence
%     set(gcf, 'Position', [100 50 width_b height_b])
% %         figure(h1)
%     ii=k;
%     curImgFrame=load(strcat(movieRawFluoPath,'/flscRaw',num2str(ii,iiformat),'.mat'));
%     curImgFrame = curImgFrame.curImgFrame;
%     imshow(curImgFrame,[]); hold on
%     set(gca,'XLim',[cropInfo1+band cropInfo3-band],'YLim',[cropInfo2+band cropInfo4-band])
% end
%% Saving imgMap and tMap in its single frames for parallel processing
disp('Saving raw mat file for parallel image overlay...')
tic
for ii=1:numFrames
    curImgFrame = imgMap(:,:,ii);
    save(strcat(movieRawFluoPath,'/flscRaw',num2str(ii,iiformat),'.mat'),'curImgFrame')
    curForceFrame = tMap(:,:,ii);
    save(strcat(movieRawForcePath,'/forceRaw',num2str(ii,iiformat),'.mat'),'curForceFrame')
end
toc
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
    try
        parpool(nWorkers)
    catch
        p= gcp;
        nWorkers = p.NumWorkers;
    end
end
% else
%     p= gcp;
%     nWorkers = p.NumWorkers;
% end
% spmd (nWorkers)
%     for k=1:floor(numFrames/nWorkers)
% h1 = figure('color','w'); % for fluorescence
% set(h1, 'Position', [100 50 width_b height_b])
% h2=figure('color','w'); % for force
% set(h2, 'Position', [100+width_b 50 width_b*1.25 height_b])
cropInfo1=cropInfo(1);
cropInfo2=cropInfo(2);
cropInfo3=cropInfo(3);
cropInfo4=cropInfo(4);
% for k=1:floor(numFrames)
% parpool(nWorkers)
parfor k=1:floor(numFrames)
    h1 = figure('color','w'); % for fluorescence
    set(h1, 'Position', [100 50 width_b height_b])
%         figure(h1)
    axIMap = axes;
%     hold off
%         ii = (k-1)*nWorkers+labindex;
    ii=k;
    curImgFrame=load(strcat(movieRawFluoPath,'/flscRaw',num2str(ii,iiformat),'.mat'));
    curImgFrame = curImgFrame.curImgFrame;
%     imshow(curImgFrame,[]); hold on
    imshow(curImgFrame,[],'Parent',axIMap); hold on
%     set(axIMap,'XLim',[cropInfo1+band cropInfo3-band],'YLim',[cropInfo2+band cropInfo4-band])
    if drawWithoutClasses==0
        for jj=1:numel(tracksNA)
            if tracksNA(jj).presence(ii)
                plot(tracksNA(jj).xCoord(ii),tracksNA(jj).yCoord(ii),'ro','MarkerSize',markerSize)
                plot(tracksNA(jj).xCoord(1:ii),tracksNA(jj).yCoord(1:ii),'r')
            end
        end
    else
        drawClassifiedTracks(allDataClass,tracksNA,ii,gca,drawWithoutClasses);
%         [htrackLine, htrackCircles] =drawClassifiedTracks(allDataClass,tracksNA,ii,gca);
    end
    set(gca,'XLim',[cropInfo1+band cropInfo3-band],'YLim',[cropInfo2+band cropInfo4-band])
    drawnow
    print('-depsc2', '-r150', strcat(movieEpsFluoPath,'/flscInten',num2str(ii,iiformat),'.eps'));
    print('-dtiff', '-r150', strcat(movieTifFluoPath,'/flscInten',num2str(ii,iiformat),'.tif'));
%     hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
%         hgsave(h2,strcat(movieFigFluoPath,'/flscInten',num2str(ii,iiformat)),'-v7.3')
    hold off
    close

%         figure(h2), 
    h2=figure('color','w'); % for force
    set(h2, 'Position', [100+width_b 50 width_b*1.25 height_b])
    subplot('Position',[0 0 0.8 1]); 
    curForceFrame=load(strcat(movieRawForcePath,'/forceRaw',num2str(ii,iiformat),'.mat'));
    curForceFrame = curForceFrame.curForceFrame;
    imshow(curForceFrame,[tmin tmax]); colormap jet; hold on;
    if drawWithoutClasses==0
        for jj=1:numel(tracksNA)
            if tracksNA(jj).presence(ii)
                plot(tracksNA(jj).xCoord(ii),tracksNA(jj).yCoord(ii),'ro','MarkerSize',markerSize)
                plot(tracksNA(jj).xCoord(1:ii),tracksNA(jj).yCoord(1:ii),'r')
            end
        end
    else
%         for kk=1:numel(htrackLine)
%             numCurLine = numel(htrackLine{kk});
%             for kkk=1:numCurLine
%                 copyobj(htrackLine{kk}{kkk},gca)
%             end
%             %drawing
% %             cellfun(@(x) copyobj(x,axTMap),htrackLine{kk})
%         end
%         for kk=1:numel(htrackCircles)
%             numCurCircle = numel(htrackCircles{kk});
%             for kkk=1:numCurCircle
%                 copyobj(htrackCircles{kk}{kkk},gca)
%             end
% %             cellfun(@(x) copyobj(x,axTMap),htrackCircles{kk})
%         end
        drawClassifiedTracks(allDataClass,tracksNA,ii,gca,drawWithoutClasses);
    end
    set(gca,'XLim',[cropInfo1+band cropInfo3-band],'YLim',[cropInfo2+band cropInfo4-band])
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
    close
end

end
% close([h1 h2])
% delete(gcp)
        