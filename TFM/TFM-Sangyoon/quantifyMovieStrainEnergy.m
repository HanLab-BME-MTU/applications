function [SE_FOV,SE_FAs,integratedForce,maxForceBlob,bioSensorAtBlob,maxForceDistFromEdge,pixelAreaForceBlob]=quantifyMovieStrainEnergy(curMD,tMax,drawFig)
%function [SE_FOV,SE_FAs,SE_CellSeg]=quantifyMovieStrainEnergy(curMD)
%quantifies from traction map the strain energy for entire field of view
%(FOV), segmented FAs, and cell segmentation when the cell segmentation
%information is there.
% Unit is in femto Joule (1e-15 J)
% Sangyoon Han March, 2016
disp(['Working on ' curMD.getFullPath '...'])
if nargin<3
    drawFig=true;
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
epsPath = [imgPath filesep 'eps'];
figPath = [imgPath filesep 'figs'];
if ~exist(figPath,'dir')
    mkdir(imgPath);
    mkdir(dataPath);
    mkdir(epsPath);
    mkdir(figPath);
end

tractionMaps=load(forceFieldProc.outFilePaths_{2});
tMap = tractionMaps.tMap; % this is currently in Pa per pixel (1pix x 1pix)
yModulus = forceFieldProc.funParams_.YoungModulus;
%% Load the displfield
iCorrectedDisplFieldProc = 3;
CorrectedDisplFieldProc=TFMPackage.processes_{iCorrectedDisplFieldProc};
if ~isempty(CorrectedDisplFieldProc)
    displMaps=load(CorrectedDisplFieldProc.outFilePaths_{2});
    dMap=displMaps.dMap; % this is currently in pix
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
maskProc = curMD.getProcess(curMD.getProcessIndex('MaskRefinementProcess'));
iSDCProc =curMD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=curMD.processes_{iSDCProc};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
end
%% Calculate strain energy for segmented FAs.
tic
minSize = 20; % in pixel
minTraction = 10; % in Pa
borderWidth=20;
nMaxForce = 0; % the number of maxForceBlobs
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
    mask = maskProc.loadChannelOutput(iChan,ii);
    if (size(mask,1)~=size(croppedTMap,1) || size(mask,2)~=size(croppedTMap,2)) && ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        % Apply subpixel-wise registration to original masks
        Ibw = padarray(mask, [maxY, maxX]);
        mask = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
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
    stats=regionprops(maskAdhesion,croppedTMap,'Area','PixelIdxList','Centroid','MaxIntensity','MeanIntensity','WeightedCentroid');
    SE_FAs.nFA(ii)=numel(stats);
    SE_FAs.areaFA(ii) = sum(maskAdhesion(:))*areaConvert; % in um2
    SE_FAs.avgFAarea(ii) = SE_FAs.areaFA(ii)/SE_FAs.nFA(ii);
    SE_FAs.avgSEperFA(ii) = SE_FAs.SE(ii)/SE_FAs.nFA(ii); % still in femto-J
    SE_FAs.SEDensity(ii)=SE_FAs.SE(ii)/SE_FAs.areaFA(ii)*1e3; % J/m2
    individualForceBlobs = arrayfun(@(x) x.MeanIntensity,stats);
    individualForceBlobMax = arrayfun(@(x) x.MaxIntensity,stats);
    individualForceBlobCenters = arrayfun(@(x) x.WeightedCentroid,stats,'UniformOutput',false);
    
    distToEdge = cellfun(@(x) min(sqrt(sum((allBdPoints- ones(size(allBdPoints,1),1)*[x(1), x(2)]).^2,2))),individualForceBlobCenters);
    integratedForce.force(ii)=sum(croppedTMap(maskAdhesion))*(pixSize_mu*1e-6)^2*1e9; % in nN
    integratedForce.avgTraction{ii}=(individualForceBlobs); % in Pa
    integratedForce.maxTraction{ii}=(individualForceBlobMax); % in Pa
    integratedForce.distToEdge{ii}=distToEdge*pixSize_mu; % in um
    integratedForce.forceBlobPixelIdxList{ii}=arrayfun(@(x) x.PixelIdxList,stats,'UniformOutput',false);
    % find an adhesion that contains max traction
    [curMaxForce,iMaxForceBlob]=max(individualForceBlobMax);
    curMaxForceBlobPixelList = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
    if ii==1
        nMaxForce=nMaxForce+1;
        maxForceBlob{nMaxForce}(ii)=curMaxForce;
        % Get the pixel list
        maxForceBlobPixelList{nMaxForce} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
    elseif ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce})) % the maxForce continues
        maxForceBlob{nMaxForce}(ii)=curMaxForce;
        maxForceBlobPixelList{nMaxForce} = union(maxForceBlobPixelList{nMaxForce},integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob});
    elseif nMaxForce>1 && ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce-1})) % the maxForce continues
        maxForceBlob{nMaxForce-1}(ii)=curMaxForce;
        maxForceBlobPixelList{nMaxForce-1} = union(maxForceBlobPixelList{nMaxForce-1},integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob});
    elseif nMaxForce>2 && ~isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce-2})) % the maxForce continues
        maxForceBlob{nMaxForce-2}(ii)=curMaxForce;
        maxForceBlobPixelList{nMaxForce-2} = union(maxForceBlobPixelList{nMaxForce-2},integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob});
    elseif isempty(intersect(curMaxForceBlobPixelList,maxForceBlobPixelList{nMaxForce})) % new maxForce appears
        nMaxForce=nMaxForce+1;
        maxForceBlobPixelList{nMaxForce} = integratedForce.forceBlobPixelIdxList{ii}{iMaxForceBlob};
    end
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

% Now I have to re-collect maxForce for entire frames
% If there is Biosensor package run, collect the ratio too
BiosensorsPackage = curMD.getPackage(curMD.getPackageIndex('BiosensorsPackage'));
ratioTiffProc = BiosensorsPackage.getProcess(11);
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
    end
    % Cell Boundary
    [B,~,nBD]  = bwboundaries(mask,'noholes');

    allBdPoints = [];
    for kk=1:nBD
        boundary = B{kk};
        allBdPoints = [allBdPoints; boundary(:,2), boundary(:,1)];
    end
    if drawFig
        if ii==1
            h1=figure; 
        else
            figure(h1), hold off
        end
        imshow(croppedTMap.*maskShrunkenBorder,[0 tMax]); colormap jet, hold on
        colorbar
        cellfun(@(x) plot(x(:,2),x(:,1), 'Color','w', 'LineWidth', 0.5),B)
    end
    for jj=1:nMaxForce
        maxForceBlob{jj}(ii) = nanmean(croppedTMap(maxForceBlobPixelList{jj}));     
        curRatioTiff=ratioTiffProc.loadChannelOutput(3,ii);
        bioSensorAtBlob{jj}(ii) = nanmean(curRatioTiff(maxForceBlobPixelList{jj}))/1000;
        curBlobMask = false(size(curRatioTiff));
%         [nRow nCol]=ind2sub(size(curRatioTiff),maxForceBlobPixelList{jj});
        curBlobMask(maxForceBlobPixelList{jj})=true;
        curBlobCentroid=regionprops(curBlobMask,'Centroid');
        distToEdge =min(sqrt(sum((allBdPoints- ones(size(allBdPoints,1),1)*[curBlobCentroid.Centroid]).^2,2)));
        maxForceDistFromEdge{jj}(ii) = distToEdge*pixSize_mu; % in um.
        adhBound = bwboundaries(curBlobMask,'noholes');    
        if drawFig
            cellfun(@(x) plot(x(:,2),x(:,1), 'Color','k', 'LineWidth', 0.5),adhBound)
        end
        if ii==1
            pixelAreaForceBlob(jj)=length(maxForceBlobPixelList{jj});
        end
    end
    if drawFig
        hgexport(h1,strcat(epsPath,'/forceMapSeg',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','eps')
        hgsave(h1,strcat(figPath,'/forceMapSeg',num2str(ii,iiformat)),'-v7.3')
        hold off
    end
end

disp('Done!')
save([dataPath filesep 'SE.mat'],'SE_FOV','SE_FAs','integratedForce','maxForceBlob','bioSensorAtBlob','maxForceDistFromEdge')

if drawFig && ~isempty(h1)
    close(h1)
end
% toc

