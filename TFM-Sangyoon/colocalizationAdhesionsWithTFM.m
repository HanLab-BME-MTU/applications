function [forceNA,forceFC,forceFA,pstruct_NAcell] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,tmax,pstruct_NAcell)
% colocalizationTFMwithFA runs colocalization between peaks in TFM maps and peaks in paxillin using movieData.
% Basically this function make grayscale of TFM and pick up peaks of the
% map and see if how many of them are colocalized with significant paxillin
% signal or neighborhood of nascent adhesions.

% input:    pathForTheMovieDataFile:    path to the movieData file (TFM
% package should be run beforehand)
%           outputPath              outputPath
%           band:                       band width for cutting border
%           (default=4)
%           pointTF:                   true if you want to trace force at
%           a certain point (default = false)
% output:   images of heatmap stored in pathForTheMovieDataFile/heatmap
%           forceNAdetec,          force at detectable NAs
%           forceNAund,             force at undetectable NAs
%           forceFA,                     force magnitude at FA
%           peaknessRatio,      portion of forces that are themselves local maxima
%           DTMS                        Deviation of Traction Magnitude from Surrounding
% Note: detectable NAs are defined by coincidence of NA position with bead
% location and force node
% Sangyoon Han April 2013

if nargin < 3
    band = 16;
    tmax=[];
    pstruct_NAcell=[];
elseif nargin <4
    tmax=[];
    pstruct_NAcell=[];
elseif nargin<5
    pstruct_NAcell=[];
end
if isempty(pstruct_NAcell)
    bFirst=true;
else
    bFirst=false;
end
bandwidth = 12; %um
%% Data Set up
% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
MD = MovieData.load(movieDataPath);
% Get whole frame number
nFrames = MD.nFrames_;
% Get TFM package
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Load the displField
iDispFieldProc = 3;
displFieldProc=TFMPackage.processes_{iDispFieldProc};
maskArray = MD.getROIMask;
% Use mask of first frame to filter displacementfield
firstMask = maskArray(:,:,1);
displFieldOriginal=displFieldProc.loadChannelOutput;
displField = filterDisplacementField(displFieldOriginal,firstMask);

% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
forceField=forceFieldProc.loadChannelOutput;

% Load the Paxillin channel

% Set up the output file path
outputFilePath = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
paxPath = [outputFilePath filesep 'pax'];
forcemapPath = [outputFilePath filesep 'fMap'];
forcetifPath = [outputFilePath filesep 'forcetifs'];
paxtifPath = [outputFilePath filesep 'paxtifs'];
epsPath = [outputFilePath filesep 'eps'];
figPath = [outputFilePath filesep 'figs'];
if ~exist(forcetifPath,'dir') || ~exist(paxtifPath,'dir') || ~exist(paxPath,'dir') || ~exist(figPath,'dir') || ~exist(epsPath,'dir') || ~exist(forcemapPath,'dir')  || ~exist(dataPath,'dir') 
    mkdir(paxPath);
    mkdir(forcetifPath);
    mkdir(paxtifPath);
    mkdir(figPath);
    mkdir(epsPath);
    mkdir(forcemapPath);
    mkdir(dataPath);
end

%Find the maximum force.
tmin = 100000;
[reg_grid1,~,~,~]=createRegGridFromDisplField(displField(1),1); %2=2 times fine interpolation

% band width for cutting border
%     band=4;
if isempty(tmax)
    tmax = 0;
    for ii = 1:nFrames
       %Load the saved body force map.
        [~,fmat, ~, ~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid1); %1:cluster size
        fnorm = (fmat(:,:,1).^2 + fmat(:,:,2).^2).^0.5;
        % Boundary cutting - I'll take care of this boundary effect later
        fnorm(end-round(band/2):end,:)=[];
        fnorm(:,end-round(band/2):end)=[];
        fnorm(1:1+round(band/2),:)=[];
        fnorm(:,1:1+round(band/2))=[];
        fnorm_vec = reshape(fnorm,[],1); 

        tmax = max(tmax,max(fnorm_vec));
        tmin = min(tmin,min(fnorm_vec));
    end
else
    for ii = 1:nFrames
       %Load the saved body force map.
        [~,fmat, ~, ~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid1); %1:cluster size
        fnorm = (fmat(:,:,1).^2 + fmat(:,:,2).^2).^0.5;
        % Boundary cutting - I'll take care of this boundary effect later
        fnorm(end-round(band/2):end,:)=[];
        fnorm(:,end-round(band/2):end)=[];
        fnorm(1:1+round(band/2),:)=[];
        fnorm(:,1:1+round(band/2))=[];
        fnorm_vec = reshape(fnorm,[],1); 

        tmin = min(tmin,min(fnorm_vec));
    end
end

[reg_grid,~,~,~]=createRegGridFromDisplField(displField(ii),4); %4=4 times fine interpolation

%     h2 = figure;
hold off
%     set(h2, 'Position', [100+imSizeX*10/9 100 imSizeX imSizeY])
hc = []; %handle for colorbar
hl = []; %handle for scale bar
iiformat = ['%.' '3' 'd'];
%     paxLevel = zeros(nFrames,1);
% SegmentationPackage = MD.getPackage(MD.getPackageIndex('SegmentationPackage'));
minSize = round((500/MD.pixelSize_)*(150/MD.pixelSize_)); %adhesion limit=1um*.5um

for ii=1:nFrames
    [~,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
    [grid_mat,if_mat,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid);

    % Boundary cutting - I'll take care of this boundary effect later
    if_mat(end-band:end,:,:)=[];
    if_mat(:,end-band:end,:)=[];
    if_mat(1:1+band,:,:)=[];
    if_mat(:,1:1+band,:)=[];
    iu_mat(end-band:end,:,:)=[];
    iu_mat(:,end-band:end,:)=[];
    iu_mat(1:1+band,:,:)=[];
    iu_mat(:,1:1+band,:)=[];
    grid_mat(end-band:end,:,:)=[];
    grid_mat(:,end-band:end,:)=[];
    grid_mat(1:1+band,:,:)=[];
    grid_mat(:,1:1+band,:)=[];

    pos = [reshape(grid_mat(:,:,1),[],1) reshape(grid_mat(:,:,2),[],1)]; %dense
    disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    force_vec = [reshape(if_mat(:,:,1),[],1) reshape(if_mat(:,:,2),[],1)]; 
    [grid_mat,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],grid_mat); %1:cluster size
    tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;

    % size of the region of interest
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1);
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2);

    % Making traction stress map, tsMap
    [XI,YI]=meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
    tsMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tnorm,XI,YI,'linear');
    % Cell Boundary Mask 
%     maskProc = SegmentationPackage.processes_{2};
%     mask = maskProc.loadChannelOutput(2,ii);
%     cropMask = mask(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    % Cell Boundary

    %%-------------Adhesion detection----------------------
    % loading paxillin image
    iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(iSDCProc)
        SDCProc=MD.processes_{iSDCProc};
        if ~SDCProc.checkChannelOutput(1)
            error(['The channel must have been corrected ! ' ...
                'Please apply stage drift correction to all needed channels before '...
                'running displacement field calclation tracking!'])
        end
        paxImage=double(SDCProc.loadChannelOutput(2,ii)); %movieData.channels_(2).loadImage(ii);
    else
        paxImage=double(MD.channels_(2).loadImage(ii)); 
    end
    paxImageCropped = paxImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);

    % Get Cell mask
    pId = double(paxImageCropped)/max(double(paxImageCropped(:)));
    % alpha = graythresh(pId);
    %estimate the intensity level to use for thresholding the image
    level1 = graythresh(pId); %Otsu
    [~, level2] = cutFirstHistMode(pId,0); %Rosin
%     alpha = 0.99*level2 + 0.01*level1;
    alpha = 0.88*level2;

    pId = double(paxImage)/max(double(paxImageCropped(:)));
    bwPI = im2bw(pId,alpha);
    bwPI2 = bwmorph(bwPI,'clean');
    bwPI3 = bwmorph(bwPI2,'erode',5);
    bwPI3 = bwmorph(bwPI3,'dilate',4);
    bwPI4 = bwmorph(bwPI3,'close',10);
    bwPI4 = bwmorph(bwPI4,'dilate',1);
    % bwPI5 = refineEdgeWithSteerableFilterGM(pId,bwPI4);
    % In case that there is still islands, pick only the largest chunk
    [labelPI,nChunk] = bwlabel(bwPI4);
    if nChunk>1
        eachArea = zeros(nChunk,1);
        for k=1:nChunk
            currBWPI = labelPI==k;
            eachArea(k) = sum(currBWPI(:));
        end
        [~,indCellArea] = max(eachArea);
        bwPI4 = labelPI == indCellArea;
    end
    cropMask = bwPI4(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    B  = bwboundaries(cropMask,'noholes');
    boundary = B{1};

    % mask for band from edge
    iMask = imcomplement(cropMask);
    distFromEdge = bwdist(iMask);
    bandwidth_pix = round(bandwidth*1000/MD.pixelSize_);
    bandMask = distFromEdge <= bandwidth_pix;

    % Get the mask for FAs
    maskAdhesion = blobSegmentThreshold(paxImageCropped,minSize,false,bandMask & cropMask);
    adhBound = bwboundaries(maskAdhesion,'noholes');
    
    % Nascent adhesion detection
    sigmaPSF_NA = getGaussianPSFsigma(MD.numAperture_, 1, MD.pixelSize_*1e-9, 601*1e-9);

    bandwidth = 3.5; %um
    bandwidth_pix = round(bandwidth*1000/MD.pixelSize_);
    bandMask = distFromEdge <= bandwidth_pix;

    naMask = bandMask & cropMask & ~maskAdhesion;
   
    % Get all nascent adhesions
    if bFirst
        pstruct_NA = pointSourceDetection(paxImageCropped, sigmaPSF_NA*1.3, 'alpha', 0.01,'Mask',naMask);

        % Showing for debug
        h0=figure; imshow(paxImageCropped,[])
        hold on;
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.25) % cell boundary
        plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
        for k = 1:length(adhBound)
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 1) %adhesion boundary
        end

        % ginput test for not-detected adhesions
        disp('Click nascent adhesions and press Enter when you are done with clicking');
        NA_man_detected = ginput;
        if ~isempty(NA_man_detected)
            plot(NA_man_detected(:,1),NA_man_detected(:,2),'ro')
            nPoints = length(pstruct_NA.A);
            old_pstruct_NA = pstruct_NA;
            pstruct_NA.x(nPoints+1:nPoints+length(NA_man_detected)) = NA_man_detected(:,1);
            pstruct_NA.y(nPoints+1:nPoints+length(NA_man_detected)) = NA_man_detected(:,2);
            % for amplitude
            for k=1:size(NA_man_detected,1)
                pstruct_NA.A(nPoints+k) = paxImageCropped(round(NA_man_detected(k,2)),round(NA_man_detected(k,1)));
            end
        end
        pstruct_NAcell{ii} = pstruct_NA;
    else
        pstruct_NA = pstruct_NAcell{ii};
    end
    
    % Showing for debug (TFM)
    h1 = figure('color','w');
    set(h1, 'Position', [100 100 (imSizeX+1)*1.25 imSizeY+1])
    subplot('Position',[0 0 0.8 1])
    imshow(tsMap,[tmin tmax]), colormap jet;
    hold on;
%     plot(pstruct.x,pstruct.y,'mo') % all peaks
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
    % unit vector plot
    [reg_grid_coarse,~,~,~]=createRegGridFromDisplField(displField(ii),1); %2=2 times fine interpolation
    [grid_mat_coarse,iu_mat_coarse,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid_coarse);
    pos_coarse = [reshape(grid_mat_coarse(:,:,1),[],1) reshape(grid_mat_coarse(:,:,2),[],1)]; %dense
    disp_vec_coarse = [reshape(iu_mat_coarse(:,:,1),[],1) reshape(iu_mat_coarse(:,:,2),[],1)]; 
    [~,if_mat_coarse,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid_coarse);
    force_vec_coarse = [reshape(if_mat_coarse(:,:,1),[],1) reshape(if_mat_coarse(:,:,2),[],1)]; 

    [~,tmat_coarse, ~, ~] = interp_vec2grid(pos_coarse+disp_vec_coarse, force_vec_coarse,[],grid_mat_coarse); %1:cluster size

    tmat_coarse(end-band/4-1:end,:,:)=[];
    tmat_coarse(:,end-band/4-1:end,:)=[];
    tmat_coarse(1:1+band/4+1,:,:)=[];
    tmat_coarse(:,1:1+band/4+1,:)=[];
    grid_mat_coarse(end-band/4-1:end,:,:)=[];
    grid_mat_coarse(:,end-band/4-1:end,:)=[];
    grid_mat_coarse(1:1+band/4+1,:,:)=[];
    grid_mat_coarse(:,1:1+band/4+1,:)=[];

    tmat_vecx = reshape(tmat_coarse(:,:,1),[],1);
    tmat_vecy = reshape(tmat_coarse(:,:,2),[],1);
    pos_vecx = reshape(grid_mat_coarse(:,:,1),[],1);
    pos_vecy = reshape(grid_mat_coarse(:,:,2),[],1);
    forceScale=0.1*max(sqrt(tmat_vecx.^2+tmat_vecy.^2));
    quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), tmat_vecx./forceScale,tmat_vecy./forceScale,0,'Color',[75/255 0/255 130/255]);

    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
    
    % Find detectable NAs from pstruct_NA
    NAs = [grid_mat(1,1,1)+pstruct_NA.x' grid_mat(1,1,2)+pstruct_NA.y'];
    deformedBeads = [displField(ii).pos(:,1)+displField(ii).vec(:,1) displField(ii).pos(:,2)+displField(ii).vec(:,2)];
    idx = KDTreeBallQuery(deformedBeads,NAs, 3);
    valid = ~cellfun(@isempty, idx);
%     valid = (displField.vec(cell2mat(idx),1).^2+displField.vec(cell2mat(idx),2).^2).^0.5;
    NAdetec = NAs(valid, :);
    NAundetec = NAs(~valid, :);

    plot(NAdetec(:,1)-grid_mat(1,1,1),NAdetec(:,2)-grid_mat(1,1,2),'ro')
    plot(NAundetec(:,1)-grid_mat(1,1,1),NAundetec(:,2)-grid_mat(1,1,2),'ro')
    forceNAdetec = zeros(length(NAdetec),1);
    neighPix = 1;
    for i=1:length(NAdetec)
        curF_NAdetec = tsMap(round(NAdetec(i,2))-grid_mat(1,1,2)-neighPix:round(NAdetec(i,2))-grid_mat(1,1,2)+neighPix,...
            round(NAdetec(i,1))-grid_mat(1,1,1)-neighPix:round(NAdetec(i,1))-grid_mat(1,1,1)+neighPix);
        forceNAdetec(i) = max(curF_NAdetec(:));
    end
    forceNAund = zeros(length(NAundetec),1);
    for i=1:length(NAundetec)
        curF_NAund = tsMap(round(NAundetec(i,2))-grid_mat(1,1,2)-neighPix:round(NAundetec(i,2))-grid_mat(1,1,2)+neighPix,...
            round(NAundetec(i,1))-grid_mat(1,1,1)-neighPix:round(NAundetec(i,1))-grid_mat(1,1,1)+neighPix);
        forceNAund(i) = max(curF_NAund(:));
    end
    forceNA = [forceNAdetec; forceNAund];
%     forceNAdetec = tsMap(round(NAdetec(:,2))-grid_mat(1,1,2),round(NAdetec(:,1))-grid_mat(1,1,1));
%     plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
%     plot(NA_man_detected(:,1),NA_man_detected(:,2),'m*')
    
    % focal contact (FC) analysis
    Adhs = regionprops(maskAdhesion,'Area','Eccentricity','PixelIdxList' );
    minFASize = round((2000/MD.pixelSize_)*(500/MD.pixelSize_)); %adhesion limit=1um*.5um

    adhIdx = arrayfun(@(x) x.Area<minFASize & x.Eccentricity<0.95, Adhs);
    FCs = Adhs(adhIdx);
    forceFC = arrayfun(@(x) max(tsMap(x.PixelIdxList)),FCs,'UniformOutput',false);
    forceFC = cell2mat(forceFC);
    FCIdx = find(adhIdx);
    for k = FCIdx'
        adhBoundary = adhBound{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'y', 'LineWidth', 0.5) %adhesion boundary
    end

    % for larger adhesions
    FAs = Adhs(~adhIdx);
    forceFA = arrayfun(@(x) max(tsMap(x.PixelIdxList)),FAs,'UniformOutput',false);
    forceFA = cell2mat(forceFA);
    FAIdx =  find(~adhIdx);
    for k = FAIdx'
        adhBoundary = adhBound{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'c', 'LineWidth', 0.5) %adhesion boundary
    end
    
    % Scale bar 2000nm
    subplot('Position',[0.8 0.1 0.1 0.8])
    axis tight
    caxis([tmin tmax]), axis off
    hc = colorbar('West');
    set(hc,'Fontsize',16)
    hold on;

    print('-depsc2', '-r300', strcat(epsPath,'/forcePeakDistribution',num2str(ii,iiformat),'.eps'));
    hgexport(h1,strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h1,strcat(figPath,'/forcePeakFig',num2str(ii,iiformat)),'-v7.3')

%     % Get magnitude of force and the curvature of the force at each (stored
%     % in pstruct_NAwithForce.fmag and pstruct_NAwithForce.fcurvature)
%     pstruct_NAwithForce = findMagCurvature(tsMap,pstruct_NA,5);
%     
%     % Plot paxillin intensity with force magnitude
%     figure, plot(pstruct_NAwithForce.A,pstruct_NAwithForce.fmag,'.')
%     xlabel('Pax Intensity (A.U.)')
%     ylabel('Stress magnitude (Pa)')
%     figure, plot(pstruct_NAwithForce.A,pstruct_NAwithForce.fcurvature,'r.')
%     % Make a histogram with pstruct_NAwithForce
%     figure, hist(pstruct_NAwithForce.fmag,250:500:3000)

%     save(strcat(dataPath,'/pstructNAForce',num2str(ii,iiformat)), 'pstruct_NAwithForce');
%     %---- for FAs, do line integral of force gradient along focal adhesion
%     distVals=-3:3; %positive =inside
%     repDist = 2.5:-1:-2.5; % making positive=outside for plotting
%     [avgFprofile,avgFprofileStd,nPts] = intensityVsDistFromEdge(tsMap,maskAdhesion,distVals);
%     figure,plot(repDist,avgFprofile)
%     [~,avgFslope,~] = regression(repDist,avgFprofile);
% %     plot(repDist,mean(avgFprofile,1),'k','LineWidth',2)
%     FAstruct.avgFprofile = avgFprofile;
%     FAstruct.avgFprofileStd = avgFprofileStd;
%     FAstruct.avgFslope = avgFslope;
%     FAstruct.Finside = avgFprofile(end);
%     FAstruct.FinsideErr = avgFprofileStd(end)/sqrt(max(nPts)); %SEM
%     save(strcat(dataPath,'/FAstruct',num2str(ii,iiformat)),'FAstruct');
    
    h2=figure;
    %Scale bar 2 um
    paxImageCropped(15:16,10:10+round(2000/MD.pixelSize_))=max(max(paxImageCropped));
    imshow(paxImageCropped,[]), hold on
    plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
    for k = 1:length(adhBound)
        adhBoundary = adhBound{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 0.5) %adhesion boundary
    end
    print('-depsc2', '-r300', strcat(epsPath,'/paxPeakDistribution',num2str(ii,iiformat),'.eps'));
    hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
   
%     close(h2)
%     close(h1)
end
end

function pstruct_NAwithForce = findMagCurvature(tsMap,pstruct_NA,neighD)
    nPoints = length(pstruct_NA.x);
    pstruct_NAwithForce = pstruct_NA;
    laplacian = [.5 1 .5; 1 -6 1; .5 1 .5];
    for jj=1:nPoints
        rowRange = round(pstruct_NA.y(jj))-neighD:round(pstruct_NA.y(jj))+neighD;
        colRange = round(pstruct_NA.x(jj))-neighD:round(pstruct_NA.x(jj))+neighD;
        pstruct_NAwithForce.fmag(jj) = max(max(tsMap(rowRange,colRange))); %force magnitude
        pstruct_NAwithForce.fcurvature(jj) = sum(sum(tsMap(round(pstruct_NA.y(jj))-1:round(pstruct_NA.y(jj))+1,round(pstruct_NA.x(jj))-1:round(pstruct_NA.x(jj))+1)...
                                                                            .* laplacian)); %force curvature
    end
end