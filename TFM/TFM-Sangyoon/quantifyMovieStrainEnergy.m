function [SE_FOV,SE_FAs,integratedForce,maxForceBlob,bioSensorAtBlob,maxForceDistFromEdge,pixelAreaForceBlob,indMaxBlob,isProtrusion]=quantifyMovieStrainEnergy(curMD,tMax,drawFig,normFactor,sRange)
%function [SE_FOV,SE_FAs,SE_CellSeg]=quantifyMovieStrainEnergy(curMD)
%quantifies from traction map the strain energy for entire field of view
%(FOV), segmented FAs, and cell segmentation when the cell segmentation
%information is there.
% Unit is in femto Joule (1e-15 J)
% Sangyoon Han March, 2016
disp(['Working on ' curMD.getFullPath '...'])
if nargin<3
    drawFig=true;
    normFactor=1;
end
if nargin<4
    normFactor=1;
end
if nargin<5
    sRange=[];
end
%% Load TFMPackage
nFrames = curMD.nFrames_;
% Get TFM package
TFMPackage = curMD.getPackage(curMD.getPackageIndex('TFMPackage'));

%% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
forceFieldPath=forceFieldProc.outFilePaths_{1};
pathFF=fileparts(forceFieldPath);
tractionImgFolder=[pathFF filesep 'tractionImgWithSegmentation'];
imgPath = [tractionImgFolder filesep 'imgs'];
dataPath = [tractionImgFolder filesep 'data'];
tifForcePath = [imgPath filesep 'tifForce'];
tifBSensorPath = [imgPath filesep 'tifBSensor'];
epsPath = [imgPath filesep 'eps'];
figPath = [imgPath filesep 'figs'];
if ~exist(tifForcePath,'dir')
    mkdir(imgPath);
    mkdir(dataPath);
    mkdir(epsPath);
    mkdir(figPath);
    mkdir(tifForcePath);
    mkdir(tifBSensorPath);
end

tractionMaps=load(forceFieldProc.outFilePaths_{2});
tMap = tractionMaps.tMap; % this is currently in Pa per pixel (1pix x 1pix)
yModulus = forceFieldProc.funParams_.YoungModulus;
%% Load the displfield
iCorrectedDisplFieldProc = 3;
CorrectedDisplFieldProc=TFMPackage.processes_{iCorrectedDisplFieldProc};
if ~isempty(CorrectedDisplFieldProc)
    try
        displMaps=load(CorrectedDisplFieldProc.outFilePaths_{2});
        dMap=displMaps.dMap; % this is currently in pix
    catch
        disp('Assuming displacement proportional to tMap because there was no iCorrectedDisplFieldProc found')
        dMap=cellfun(@(x) x/(yModulus),tMap,'UniformOutput',false);
    end
else
    disp('Assuming displacement proportional to tMap because there was no iCorrectedDisplFieldProc found')
    dMap=cellfun(@(x) x/(yModulus),tMap,'UniformOutput',false);
end
%% Calculate strain energy for FOV
% gridSpacing=1;
pixSize_mu=curMD.pixelSize_*1e-3;
% factorConvert=gridSpacing^2*pixSize_mu^3/10^6;
% factorConvert=gridSpacing^2*(pixSize_mu*1e-6)^3;
areaConvert=pixSize_mu^2;
SE_FOV=struct('SE',zeros(nFrames,1),'area',[],'SEDensity',[]);
SE_FAs=struct('SE',[],'nFA',[],'areaFA',[],'SEDensity',[]);
iiformat=['%.' '3' 'd'];
%% Get the cell boundary
% Cell Boundary Mask 
iChan=2;
iBeadChan=1;
iSDCProc =curMD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=curMD.processes_{iSDCProc};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
end
if ~isempty(normFactor)
    BiosensorsPackage = curMD.getPackage(curMD.getPackageIndex('BiosensorsPackage'));
    ratioTiffProc = BiosensorsPackage.getProcess(10);
    % maskProc = BiosensorsPackage.getProcess(5);
end
if ~isempty(normFactor)
    maskProc = curMD.getProcess(curMD.getProcessIndex('MaskIntersectionProcess'));
else
    maskProc = curMD.getProcess(curMD.getProcessIndex('MaskRefinementProcess'));
end
%% Calculate strain energy for segmented FAs.
tic
minSize = 20; % in pixel
minTraction = 10; % in Pa
borderWidth=40;
nMaxForce = 0; % the number of maxForceBlobs
nTopBlobs=20; % the number of top maxForceBlobs in single frames

for ii=1:nFrames
    % segment
    maskNoNaN = ~isnan(tMap{ii});
    % To crop, you can calculate the top, bottom, left, and right columns of the mask's "true" area by using any(). 
    anyCol=any(maskNoNaN);
    anyRow=any(maskNoNaN,2);
    row1 = find(anyRow,1);
    row2 = find(anyRow,1,'last');
    col1 = find(anyCol,1);
    col2 = find(anyCol,1,'last');
    
    croppedTMap = tMap{ii}(row1:row2, col1:col2);
    croppedDMap = dMap{ii}(row1:row2, col1:col2);

    % Shift cell mask according to SDC
    mask = maskProc.loadChannelOutput(2,ii);
    if (size(mask,1)~=size(croppedTMap,1) || size(mask,2)~=size(croppedTMap,2)) && ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        % Apply subpixel-wise registration to original masks
        Ibw = padarray(mask, [maxY, maxX]);
        mask = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
        if (size(mask,1)>size(croppedTMap,1) || size(mask,2)>size(croppedTMap,2))
            mask = mask(1:size(croppedTMap,1),1:size(croppedTMap,2));
        end
    end
    % Cell Boundary
    [B,~,nBD]  = bwboundaries(mask,'noholes');

    allBdPoints = [];
    for kk=1:nBD
        boundary = B{kk};
        allBdPoints = [allBdPoints; boundary(:,2), boundary(:,1)];
    end
    
    maskShrunkenBorder = true(size(croppedTMap));
    maskShrunkenBorder = bwmorph(maskShrunkenBorder,'erode',borderWidth);

    maskAdhesion = blobSegmentThreshold(croppedTMap(borderWidth+1:end-borderWidth,borderWidth+1:end-borderWidth),minSize,0,mask(borderWidth+1:end-borderWidth,borderWidth+1:end-borderWidth));
    maskAdhesion = bwmorph(maskAdhesion,'dilate',2);
    maskAdhesion = padarray(maskAdhesion,[borderWidth borderWidth]);
    maskHighTraction=croppedTMap>minTraction;
    maskAdhesion = maskAdhesion & maskHighTraction;
    
    SE_FOV.SE(ii)=1/2*sum(sum(croppedDMap.*croppedTMap))*(pixSize_mu*1e-6)^3; % this is in Newton*m. %*factorConvert;
    SE_FOV.SE(ii)=SE_FOV.SE(ii)*1e15; % this is now in femto-Joule
    SE_FOV.area(ii)=sum(sum(~isnan(croppedTMap)))*areaConvert; % this is in um2
    SE_FOV.SEDensity(ii)=SE_FOV.SE(ii)/SE_FOV.area(ii)*1e3; % J/m2
    SE_FAs.SE(ii)=1/2*sum(croppedDMap(maskAdhesion).*croppedTMap(maskAdhesion))*(pixSize_mu*1e-6)^3; % this is in Newton*m.
    SE_FAs.SE(ii)=SE_FAs.SE(ii)*1e15; % this is now in femto-Joule
    stats=regionprops(maskAdhesion,croppedTMap,'Area','PixelIdxList','Centroid','MinIntensity','MaxIntensity','MeanIntensity','WeightedCentroid');
    SE_FAs.nFA(ii)=numel(stats);
    SE_FAs.areaFA(ii) = sum(maskAdhesion(:))*areaConvert; % in um2
    SE_FAs.avgFAarea(ii) = SE_FAs.areaFA(ii)/SE_FAs.nFA(ii);
    SE_FAs.avgSEperFA(ii) = SE_FAs.SE(ii)/SE_FAs.nFA(ii); % still in femto-J
    SE_FAs.SEDensity(ii)=SE_FAs.SE(ii)/SE_FAs.areaFA(ii)*1e3; % J/m2
    individualForceBlobs = arrayfun(@(x) x.MeanIntensity,stats);
    individualForceBlobMax = arrayfun(@(x) x.MaxIntensity,stats);
%     individualForceBlobMin = arrayfun(@(x) x.MinIntensity,stats);
    individualForceBlobCenters = arrayfun(@(x) x.WeightedCentroid,stats,'UniformOutput',false);
    individualForceBlobAreas = arrayfun(@(x) x.Area,stats);
    
    distToEdge = cellfun(@(x) min(sqrt(sum((allBdPoints- ones(size(allBdPoints,1),1)*[x(1), x(2)]).^2,2))),individualForceBlobCenters);
    integratedForce.force(ii)=sum(croppedTMap(maskAdhesion))*(pixSize_mu*1e-6)^2*1e9; % in nN
    integratedForce.avgTraction{ii}=(individualForceBlobs); % in Pa
    integratedForce.maxTraction{ii}=(individualForceBlobMax); % in Pa
    integratedForce.distToEdge{ii}=distToEdge*pixSize_mu; % in um
    integratedForce.forceBlobPixelIdxList{ii}=arrayfun(@(x) x.PixelIdxList,stats,'UniformOutput',false);
    % find an adhesion that contains top three max traction
%     [individualForceBlobMaxSorted, topIDs]=sort(individualForceBlobMax,'descend');
%     nTopBlobs = min(nTopBlobs,length(individualForceBlobMax));
%     iMaxForceBlobTop = topIDs(1:nTopBlobs);
%     curMaxForceTop = individualForceBlobMaxSorted(1:nTopBlobs);
    % Switching to area-based criteria
    [individualForceBlobAreaSorted, topIDs]=sort(individualForceBlobAreas,'descend');
    nTopBlobs = min(nTopBlobs,length(individualForceBlobAreas));
    iMaxForceBlobTop = topIDs(1:nTopBlobs);
    curMaxForceTop = individualForceBlobAreaSorted(1:nTopBlobs);
    newPixels=integratedForce.forceBlobPixelIdxList{ii}(topIDs(1:nTopBlobs));
    newMaskAdhesion = false(size(maskAdhesion));
    for qq=1:numel(newPixels)
        newMaskAdhesion(newPixels{qq})=true;
    end
    
    % Create interior mask: 5 um smaller than cell mask
%     interiorMask=bwmorph(mask,'erode',5/pixSize_mu);
    distMask = bwdist(~mask);
    interiorMask = distMask>=5/pixSize_mu;
    
    peripheralMask=mask & ~interiorMask;
    nForceBlobInterior=0;
    nForceBlobPeriphery=0;
%     [curMaxForce,iMaxForceBlob]=max(individualForceBlobMax);
    for pp=1:nTopBlobs
        iMaxForceBlob = iMaxForceBlobTop(pp);
%         curMaxForce = curMaxForceTop(pp);
        curMaxForce = integratedForce.maxTraction{ii}(iMaxForceBlob);
        curMaxForceBlobPixelList = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
        if ii==1
            nMaxForce=nMaxForce+1;
            maxForceBlob{nMaxForce}(ii)=curMaxForce;
            % Get the pixel list
            maxForceBlobPixelList{nMaxForce} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
    %     elseif ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce}{min(ii-1,numel(maxForceBlobPixelList{nMaxForce}))})) % the maxForce continues
    %         maxForceBlob{nMaxForce}(ii)=curMaxForce;
    % %         maxForceBlobPixelList{nMaxForce}{ii} = union(maxForceBlobPixelList{nMaxForce},integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob});
    %         maxForceBlobPixelList{nMaxForce}{ii} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
        else %if nMaxForce>=nTopBlobs 
            iPrev=nMaxForce-1;
            matchedBlobFound=false;
            while ~matchedBlobFound && iPrev>=0 % Run through all the force blobs and seek matching blobs with current max blobs
    %             if ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce-iPrev}{min(ii-1,numel(maxForceBlobPixelList{nMaxForce-iPrev}))})) % the maxForce continues
                if ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce-iPrev})) % the maxForce continues
                    maxForceBlob{nMaxForce-iPrev}(ii)=curMaxForce;
                    maxForceBlobPixelList{nMaxForce-iPrev} = union(maxForceBlobPixelList{nMaxForce-iPrev},integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob});
                    matchedBlobFound= true;
                end
                %                 maxForceBlobPixelList{nMaxForce-iPrev}{ii} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
                iPrev=iPrev-1;
            end
    %     elseif nMaxForce>2 && ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce-2}{min(ii-1,numel(maxForceBlobPixelList{nMaxForce-2}))})) % the maxForce continues
    %         maxForceBlob{nMaxForce-2}(ii)=curMaxForce;
    % %         maxForceBlobPixelList{nMaxForce-2}{ii} = union(maxForceBlobPixelList{nMaxForce-2},integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob});
    %         maxForceBlobPixelList{nMaxForce-2}{ii} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
    %     elseif isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce}{min(ii-1,numel(maxForceBlobPixelList{nMaxForce}))})) % new maxForce appears
            if ~matchedBlobFound %&& isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce})) % new maxForce appears
                nMaxForce=nMaxForce+1;
                maxForceBlob{nMaxForce}(ii)=curMaxForce;
                maxForceBlobPixelList{nMaxForce} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
            end
        end
        % Separating interior force blobs and peripheral ones
        if sum(interiorMask(curMaxForceBlobPixelList))>sum(peripheralMask(curMaxForceBlobPixelList))
            nForceBlobInterior=nForceBlobInterior+1;
            % For debug,
            compImg(:,:,1)=interiorMask;
            newMaskAdhesion = false(size(maskAdhesion));
            newMaskAdhesion(curMaxForceBlobPixelList)=true;
            compImg(:,:,2)=newMaskAdhesion;
            compImg(:,:,3)=peripheralMask;
            figure, imshow(double(compImg))
        else
            nForceBlobPeriphery=nForceBlobPeriphery+1;
        end
    end
    % See if there is overlap between interiorMask and maskAdhesion - for
    % this I have to use AND operator
    
    
    for k=1:SE_FAs.nFA(ii)
        SE_FAs.individualSE(ii).SE(k)=1/2*sum(croppedDMap(stats(k).PixelIdxList).*croppedTMap(stats(k).PixelIdxList))*(pixSize_mu*1e-6)^3*1e15;
        SE_FAs.individualSE(ii).area(k)=stats(k).Area*areaConvert;
        SE_FAs.individualSE(ii).SED(k)=SE_FAs.individualSE(ii).SE(k)/SE_FAs.individualSE(ii).area(k)*1e3;% J/m2.
    end
    [maxT]=max(croppedTMap(:));
    SE_FAs.maxT(ii) = maxT;
    [SE_FAs.maxSEFA(ii),maxSE_FA]=max(SE_FAs.individualSE(ii).SE);
    [SE_FAs.maxSED_FA(ii),maxSED_FA]=max(SE_FAs.individualSE(ii).SED);
    
%     if nargin<2
%         tMax = max(croppedTMap(:).*maskShrunkenBorder(:));
%     end
%     if ii==1
%         h1=figure; 
%     else
%         figure(h1), hold off
%     end
%     imshow(croppedTMap.*maskShrunkenBorder,[0 tMax]); colormap jet, hold on
%     colorbar
%     adhBound = bwboundaries(maskAdhesion,'noholes');    
%     cellfun(@(x) plot(x(:,2),x(:,1), 'Color','k', 'LineWidth', 0.5),adhBound)
%     cellfun(@(x) plot(x(:,2),x(:,1), 'Color','w', 'LineWidth', 0.5),B)
%     plot(adhBound{maxSE_FA}(:,2),adhBound{maxSE_FA}(:,1), 'Color','w', 'LineWidth', 2)
%     plot(adhBound{maxSED_FA}(:,2),adhBound{maxSED_FA}(:,1), 'Color','k', 'LineWidth', 1)
%     text(adhBound{maxSE_FA}(1,2),adhBound{maxSE_FA}(1,1)-60,{'Segmented total'; ['strain energy: ' num2str(SE_FAs.maxSEFA(ii)) ' fJ']},'Color','w')
%     text(adhBound{maxSED_FA}(1,2),adhBound{maxSED_FA}(1,1)+60,{'Strain energy'; ['density: ' num2str(SE_FAs.maxSED_FA(ii)) ' J/m^2']},'Color','w')
%     text(40,80,{'Total strain energy'; ['in all segmentations: ' num2str(SE_FAs.SE(ii)) ' fJ']},'Color','w')
%     text(40,row2-80,{'Strain energy density'; ['in all segmentations: ' num2str(SE_FAs.SEDensity(ii)) ' J/m^2']},'Color','w')
% 
%     hgexport(h1,strcat(epsPath,'/forceMapSeg',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','eps')
%     hgsave(h1,strcat(figPath,'/forceMapSeg',num2str(ii,iiformat)),'-v7.3')
%     hold off
%     h2=figure; plot(distToEdge*pixSize_mu,individualForceBlobs,'ro'), 
%     hgexport(h2,strcat(epsPath,'/forceMapSeg',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','eps')
%     hgsave(h2,strcat(figPath,'/forceMapSeg',num2str(ii,iiformat)),'-v7.3')
%     close
%     h3=figure, histogram(individualForceBlobs), close
end
% Re-sort out maxForceBlobPixelList
redundantBlobs=[];
for jj=2:nMaxForce
    curMaxForceBlobPixelList = maxForceBlobPixelList{jj};
    for rr=1:jj-1
        if ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{rr})) % the maxForce continues
            maxForceBlobPixelList{rr} = union(maxForceBlobPixelList{rr},curMaxForceBlobPixelList);
            redundantBlobs= [redundantBlobs jj];
            break
        end
    end
end
maxForceBlobPixelList(redundantBlobs)=[];
maxForceBlob(redundantBlobs)=[];
nMaxForce = nMaxForce-length(redundantBlobs);
% Now I have to re-collect maxForce for entire frames
% If there is Biosensor package run, collect the ratio too

for ii=1:nFrames
    % segment
    maskNoNaN = ~isnan(tMap{ii});
    % To crop, you can calculate the top, bottom, left, and right columns of the mask's "true" area by using any(). 
    anyCol=any(maskNoNaN);
    anyRow=any(maskNoNaN,2);
    row1 = find(anyRow,1);
    row2 = find(anyRow,1,'last');
    col1 = find(anyCol,1);
    col2 = find(anyCol,1,'last');
    
    croppedTMap = tMap{ii}(row1:row2, col1:col2);

    % Shift cell mask according to SDC
    mask = maskProc.loadChannelOutput(iChan,ii);
    if (size(mask,1)~=size(croppedTMap,1) || size(mask,2)~=size(croppedTMap,2)) && ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        % Apply subpixel-wise registration to original masks
        Ibw = padarray(mask, [maxY, maxX]);
        mask = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
        if (size(mask,1)>size(croppedTMap,1) || size(mask,2)>size(croppedTMap,2))
            mask = mask(1:size(croppedTMap,1),1:size(croppedTMap,2));
        end
    end
    % Cell Boundary
    [B,~,nBD]  = bwboundaries(mask,'noholes');

    allBdPoints = [];
    for kk=1:nBD
        boundary = B{kk};
        allBdPoints = [allBdPoints; boundary(:,2), boundary(:,1)];
    end
    if ~isempty(normFactor)
        curRatioTiffNormalized=ratioTiffProc.loadChannelOutput(3,ii)/normFactor;
    end
    
    if drawFig
        if ii==1
            h1=figure; 
            h1.Position=[100 400 400 320];
            if ~isempty(normFactor)
                h2=figure; 
                h2.Position=[700 400 400 320];
            end
%         else
%             figure(h1), hold off
        end
        figure(h1), hold off
        imshow(croppedTMap.*maskShrunkenBorder,[0 tMax]); colormap jet, hold on
        colorbar
        cellfun(@(x) plot(x(:,2),x(:,1), 'Color','w', 'LineWidth', 0.5),B)
        if ~isempty(normFactor)
            figure(h2), hold off
            curRatioTiffFiltered= filterGauss2D(curRatioTiffNormalized,1);
            imshow(curRatioTiffFiltered,sRange); colormap jet, hold on
            colorbar
        end
%         cellfun(@(x) plot(x(:,2),x(:,1), 'Color','w', 'LineWidth', 0.5),B)

    end
    for jj=1:nMaxForce
%         if ii>numel(maxForceBlobPixelList{jj}) || isempty(maxForceBlobPixelList{jj}{ii})
            % back-calculate previous force blobs
            % find out frame number that has the pixels
            % No, take the maximum one
            firstPresencePixel = maxForceBlobPixelList{jj};
%             nfirstPresencePixel=find(~cellfun(@isempty,maxForceBlobPixelList{jj}),1);
%             firstPresencePixel = maxForceBlobPixelList{jj}{nfirstPresencePixel};
            % compare with cell mask
            curBlobMask = false(size(croppedTMap));
    %         [nRow nCol]=ind2sub(size(curRatioTiff),maxForceBlobPixelList{jj});
            curBlobMask(firstPresencePixel)=true;
            curBlobMask = curBlobMask & mask;
            % fill hole
            curBlobMask = imfill(curBlobMask,'holes');
            curBlobCentroid=regionprops(curBlobMask,'Centroid','Area','PixelIdxList');
            if numel(curBlobCentroid)>1
                [~,iMaxArea]=max(arrayfun(@(x) x.Area,curBlobCentroid));
                curBlobCentroid = curBlobCentroid(iMaxArea);
            end
            if isempty(curBlobCentroid)
                maxForceBlob{jj}(ii) = NaN;     
                bioSensorAtBlob{jj}(ii) = NaN;
                maxForceDistFromEdge{jj}(ii) = NaN; % in um.
                pixelAreaForceBlob(jj,ii)=NaN;
                continue
            else
                maxForceBlobPixelListNew{jj}{ii}=curBlobCentroid.PixelIdxList;
            end
%         else
%             curBlobMask = false(size(curRatioTiffNormalized));
%     %         [nRow nCol]=ind2sub(size(curRatioTiff),maxForceBlobPixelList{jj});
%             curBlobMask(maxForceBlobPixelList{jj}{ii})=true;
%             curBlobMask = curBlobMask & mask;
%             curBlobCentroid=regionprops(curBlobMask,'Centroid','Area','PixelIdxList');
%             if numel(curBlobCentroid)>1
%                 [~,iMaxArea]=max(arrayfun(@(x) x.Area,curBlobCentroid));
%                 curBlobCentroid = curBlobCentroid(iMaxArea);
%             end
%         end
        maxForceBlob{jj}(ii) = nanmean(croppedTMap(maxForceBlobPixelListNew{jj}{ii}));     
        if ~isempty(normFactor)
            bioSensorAtBlob{jj}(ii) = nanmean(curRatioTiffNormalized(maxForceBlobPixelListNew{jj}{ii}));
        else
            bioSensorAtBlob{jj}(ii) = NaN;
        end
        [distToEdge,indMinBdPoint] =min(sqrt(sum((allBdPoints- ones(size(allBdPoints,1),1)*[curBlobCentroid.Centroid]).^2,2)));
        
        closestBdPoint{jj}{ii} = allBdPoints(indMinBdPoint,:); % this is lab frame of reference. (not relative to adhesion position)
        centroidBlob{jj}{ii}=[curBlobCentroid.Centroid];
        
        maxForceDistFromEdge{jj}(ii) = distToEdge*pixSize_mu; % in um.
        pixelAreaForceBlob(jj,ii)=curBlobCentroid.Area;
    end
    % is this blob this time max?
    [~,curIndMaxBlob]=nanmax(cellfun(@(x) x(ii),maxForceBlob));
    indMaxBlob(ii) = curIndMaxBlob;
    
    if drawFig
        idxNoNan = ~cellfun(@(x) isnan(x(ii)), maxForceBlob);
        for jj=find(idxNoNan)
            curBlobMask = false(size(croppedTMap));
            curBlobMask(maxForceBlobPixelListNew{jj}{ii})=true;
            adhBound = bwboundaries(curBlobMask,'noholes');    
            figure(h1), hold on
            if ii==1 %Just to make sure
                colormap jet
                colorbar                
            end
            if jj==indMaxBlob(ii)
                cellfun(@(x) plot(x(:,2),x(:,1), 'Color','m', 'LineWidth', 2),adhBound)
            else
                cellfun(@(x) plot(x(:,2),x(:,1), 'Color','m', 'LineWidth', 0.5),adhBound)
            end
            plot(centroidBlob{jj}{ii}(1),centroidBlob{jj}{ii}(2),'mo')
            line([centroidBlob{jj}{ii}(1) closestBdPoint{jj}{ii}(1)],[centroidBlob{jj}{ii}(2) closestBdPoint{jj}{ii}(2)],'Color','y')
            text(centroidBlob{jj}{ii}(1),centroidBlob{jj}{ii}(2)+30,['ID: ' num2str(jj)],'Color','w')
            text((centroidBlob{jj}{ii}(1)+closestBdPoint{jj}{ii}(1))/2,(centroidBlob{jj}{ii}(2)+closestBdPoint{jj}{ii}(2))/2+12,...
                ['D: ' num2str(maxForceDistFromEdge{jj}(ii))],'Color','w')
            if ~isempty(normFactor)
                figure(h2), hold on
                if jj==indMaxBlob(ii)
                    cellfun(@(x) plot(x(:,2),x(:,1), 'Color','m', 'LineWidth', 2),adhBound)
                else
                    cellfun(@(x) plot(x(:,2),x(:,1), 'Color','m', 'LineWidth', 0.5),adhBound)
                end
                plot(centroidBlob{jj}{ii}(1),centroidBlob{jj}{ii}(2),'mo')
                line([centroidBlob{jj}{ii}(1) closestBdPoint{jj}{ii}(1)],[centroidBlob{jj}{ii}(2) closestBdPoint{jj}{ii}(2)],'Color','y')
                text(centroidBlob{jj}{ii}(1),centroidBlob{jj}{ii}(2)+30,['ID: ' num2str(jj)],'Color','w')
                text((centroidBlob{jj}{ii}(1)+closestBdPoint{jj}{ii}(1))/2,(centroidBlob{jj}{ii}(2)+closestBdPoint{jj}{ii}(2))/2+12,...
                    ['D: ' num2str(maxForceDistFromEdge{jj}(ii))],'Color','w')
            end
        end
%         hgexport(h1,strcat(epsPath,'/forceMapSeg',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','eps')
%         print(h1, strcat(tifForcePath,'/forceMapSeg',num2str(ii,iiformat),'.tif'),'-dtiff', '-r150');
        figure(h1)
        export_fig(strcat(epsPath,'/forceMapSeg',num2str(ii,iiformat),'.eps'));
        export_fig(strcat(tifForcePath,'/forceMapSeg',num2str(ii,iiformat),'.tif'));
%         hgsave(h1,strcat(figPath,'/forceMapSeg',num2str(ii,iiformat)),'-v7.3')
%         hgexport(h2,strcat(epsPath,'/bSensorMapSeg',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','eps')
%         print(h2, strcat(tifBSensorPath,'/bSensorMapSeg',num2str(ii,iiformat),'.tif'),'-dtiff', '-r150');
%         hgsave(h2,strcat(figPath,'/bSensorMapSeg',num2str(ii,iiformat)),'-v7.3')
        if ~isempty(normFactor)
            figure(h2)
            export_fig(strcat(epsPath,'/bSensorMapSeg',num2str(ii,iiformat),'.eps'));
            export_fig(strcat(tifBSensorPath,'/bSensorMapSeg',num2str(ii,iiformat),'.tif'));
        end
        hold off
    end
end
disp('Protrusion / retraction analysis...')
deltaT=curMD.timeInterval_;
cropMaskStack = false(size(mask,1),size(mask,2),nFrames);
%% Matching with segmented adhesions - do only when nFrames>2
if nFrames>2
    for ii=1:nFrames
        % Cell Boundary Mask 
        cropMaskStack(:,:,ii) = maskProc.loadChannelOutput(iChan,ii);
    end
    isProtrusion = false(1,nMaxForce);
    for jj=1:nMaxForce
        % Protrusion or retraction?
        nFinalFrameBlob = length(closestBdPoint{jj});
        presIdx = true(1,nFinalFrameBlob);
        % get the instantaneous velocity
        % get the distance first
        curTrackVelLength=sum(presIdx)-1;
        distTrajec=zeros(curTrackVelLength,1);
        presIdxSeq = find(presIdx);
        for kk=1:curTrackVelLength
            real_kk = presIdxSeq(kk);
            if ~isempty(closestBdPoint{jj}{real_kk+1}) && ~isempty(closestBdPoint{jj}{real_kk})
                distTrajec(kk) = sqrt(sum((closestBdPoint{jj}{real_kk+1} - ...
                    closestBdPoint{jj}{real_kk}).^2,2));
                lastPointIntX = round(closestBdPoint{jj}{real_kk+1}(1));
                lastPointIntY = round(closestBdPoint{jj}{real_kk+1}(2));
                if cropMaskStack(lastPointIntY,lastPointIntX,real_kk) %if the last point is in the first mask, it is inward
                    distTrajec(kk) = -distTrajec(kk);
                end
            else
                distTrajec(kk) = NaN;
                lastPointIntX = NaN;
                lastPointIntY = NaN;
            end
        end
        if any(distTrajec~=0)
            [Protrusion,Retraction] = getPersistenceTime(distTrajec,deltaT);%,'plotYes',true)
            if any(isnan(Retraction.persTime)) || sum(Protrusion.persTime) - sum(Retraction.persTime)>0 % this is protrusion for this track
                isProtrusion(jj) = true;
            else
                isProtrusion(jj) = false;
            end
        end
    end
    save([dataPath filesep 'SE.mat'],'SE_FOV','SE_FAs','integratedForce','maxForceBlob','bioSensorAtBlob','maxForceDistFromEdge','pixelAreaForceBlob','indMaxBlob','isProtrusion')
else
    save([dataPath filesep 'SE.mat'],'SE_FOV','SE_FAs','integratedForce','maxForceBlob','bioSensorAtBlob','maxForceDistFromEdge','pixelAreaForceBlob','indMaxBlob')
end
disp('Done!')

if drawFig && ~isempty(h1)
    close(h1)
    if ~isempty(normFactor)
        close(h2)
    end
end
% toc

