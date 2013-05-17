function [pstruct_NAcell] = colocalizationAdhesionsWithTFM( pathForTheMovieDataFile,outputPath,band,pointTF,tmax,pstruct_NAcell )
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
% Sangyoon Han April 2013

if nargin < 3
    band = 4;
    pointTF = false;
    tmax=[];
    pstruct_NAcell=[];
elseif nargin <4
    pointTF = false;
    tmax=[];
    pstruct_NAcell=[];
elseif nargin<5
    tmax=[];
    pstruct_NAcell=[];
elseif nargin<6
    pstruct_NAcell=[];
end
if isempty(pstruct_NAcell)
    bFirst=true;
else
    bFirst=false;
end
bandwidth = 6; %um
%% Data Set up
% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
movieData = MovieData.load(movieDataPath);
% Get whole frame number
nFrames = movieData.nFrames_;
% Get TFM package
TFMPackage = movieData.getPackage(movieData.getPackageIndex('TFMPackage'));
% Load the displField
iDispFieldProc = 3;
displFieldProc=TFMPackage.processes_{iDispFieldProc};
maskArray = movieData.getROIMask;
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
[reg_grid1,~,~,~]=createRegGridFromDisplField(displField,1); %2=2 times fine interpolation

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

[reg_grid,~,~,~]=createRegGridFromDisplField(displField,4); %4=4 times fine interpolation

%     h2 = figure;
hold off
%     set(h2, 'Position', [100+imSizeX*10/9 100 imSizeX imSizeY])
hc = []; %handle for colorbar
hl = []; %handle for scale bar
iiformat = ['%.' '3' 'd'];
%     paxLevel = zeros(nFrames,1);
SegmentationPackage = movieData.getPackage(movieData.getPackageIndex('SegmentationPackage'));
minSize = round((1000/movieData.pixelSize_)*(200/movieData.pixelSize_)); %adhesion limit=1um*.5um

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
    tsMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tnorm,XI,YI,'cubic');
    % Cell Boundary Mask 
    maskProc = SegmentationPackage.processes_{2};
    mask = maskProc.loadChannelOutput(2,ii);
    cropMask = mask(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    % Cell Boundary
    B  = bwboundaries(cropMask,'noholes');
    boundary = B{1};

    %%-------------Adhesion detection----------------------
    % loading paxillin image
    iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(iSDCProc)
        SDCProc=movieData.processes_{iSDCProc};
        if ~SDCProc.checkChannelOutput(1)
            error(['The channel must have been corrected ! ' ...
                'Please apply stage drift correction to all needed channels before '...
                'running displacement field calclation tracking!'])
        end
        paxImage=double(SDCProc.loadChannelOutput(2,ii)); %movieData.channels_(2).loadImage(ii);
    end
    paxImageCropped = paxImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);

    % Get the mask for FAs
    maskAdhesion = blobSegmentThreshold(paxImageCropped,minSize,false,cropMask);
    adhBound = bwboundaries(maskAdhesion,'noholes');
    
    % Nascent adhesion detection
    sigmaPSF_NA = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, 601*1e-9);

    % mask for band from edge
    iMask = imcomplement(cropMask);
    distFromEdge = bwdist(iMask);
    bandwidth_pix = round(bandwidth*1000/movieData.pixelSize_);
    bandMask = distFromEdge <= bandwidth_pix;
        
    naMask = bandMask & cropMask & ~maskAdhesion;
   
    % Get all nascent adhesions
    if bFirst
        pstruct_NA = pointSourceDetection(paxImageCropped, sigmaPSF_NA*2, 'alpha', 0.01,'Mask',naMask);

        % Showing for debug
        figure; imshow(paxImageCropped,[])
        hold on;
    %     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
        plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
        for k = 1:length(adhBound)
            adhBoundary = adhBound{k};
            plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 1) %adhesion boundary
        end

        % ginput test for not-detected adhesions
        disp('Click nascent adhesions and press Enter when you are done with clicking');
        NA_man_detected = ginput;
        plot(NA_man_detected(:,1),NA_man_detected(:,2),'m*')
        nPoints = length(pstruct_NA.A);
        old_pstruct_NA = pstruct_NA;
        pstruct_NA.x(nPoints+1:nPoints+length(NA_man_detected)) = NA_man_detected(:,1);
        pstruct_NA.y(nPoints+1:nPoints+length(NA_man_detected)) = NA_man_detected(:,2);
        % for amplitude
        for k=1:length(NA_man_detected)
            pstruct_NA.A(nPoints+k) = paxImageCropped(round(NA_man_detected(k,2)),round(NA_man_detected(k,1)));
        end
        pstruct_NAcell{ii} = pstruct_NA;
    else
        pstruct_NA = pstruct_NAcell{ii};
    end
    
    % Showing for debug (TFM)
    h1 = figure('color','w');
    set(h1, 'Position', [100 100 imSizeX*1.25 imSizeY])
    subplot('Position',[0 0 0.8 1])
    imshow(tsMap,[tmin tmax]), colormap jet;
    hold on;
%     plot(pstruct.x,pstruct.y,'mo') % all peaks
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
    % unit vector plot
    [reg_grid_coarse,~,~,~]=createRegGridFromDisplField(displField,1); %2=2 times fine interpolation
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
    forceScale=0.25*sqrt(tmat_vecx.^2+tmat_vecy.^2);
    quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), tmat_vecx./forceScale,tmat_vecy./forceScale,0,'Color',[75/255 0/255 130/255]);

    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
    plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
%     plot(NA_man_detected(:,1),NA_man_detected(:,2),'m*')
    for k = 1:length(adhBound)
        adhBoundary = adhBound{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'k', 'LineWidth', 1) %adhesion boundary
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

    % Get magnitude of force and the curvature of the force at each (stored
    % in pstruct_NAwithForce.fmag and pstruct_NAwithForce.fcurvature)
    pstruct_NAwithForce = findMagCurvature(tsMap,pstruct_NA,5);
    
%     % Plot paxillin intensity with force magnitude
%     figure, plot(pstruct_NAwithForce.A,pstruct_NAwithForce.fmag,'.')
%     xlabel('Pax Intensity (A.U.)')
%     ylabel('Stress magnitude (Pa)')
%     figure, plot(pstruct_NAwithForce.A,pstruct_NAwithForce.fcurvature,'r.')
%     % Make a histogram with pstruct_NAwithForce
%     figure, hist(pstruct_NAwithForce.fmag,250:500:3000)

    save(strcat(dataPath,'/pstructNAForce',num2str(ii,iiformat)), 'pstruct_NAwithForce');
    %---- for FAs, do line integral of force gradient along focal adhesion
    distVals=-3:3; %positive =inside
    repDist = 2.5:-1:-2.5; % making positive=outside for plotting
    [avgFprofile,avgFprofileStd,nPts] = intensityVsDistFromEdge(tsMap,maskAdhesion,distVals);
    figure,plot(repDist,avgFprofile)
    [~,avgFslope,~] = regression(repDist,avgFprofile);
%     plot(repDist,mean(avgFprofile,1),'k','LineWidth',2)
    FAstruct.avgFprofile = avgFprofile;
    FAstruct.avgFprofileStd = avgFprofileStd;
    FAstruct.avgFslope = avgFslope;
    FAstruct.Finside = avgFprofile(end);
    FAstruct.FinsideErr = avgFprofileStd(end)/sqrt(max(nPts)); %SEM
    save(strcat(dataPath,'/FAstruct',num2str(ii,iiformat)),'FAstruct');
    
    h2=figure;
    %Scale bar 2 um
    paxImageCropped(15:16,10:10+round(2000/movieData.pixelSize_))=max(max(paxImageCropped));
    imshow(paxImageCropped,[]), hold on
    plot(pstruct_NA.x,pstruct_NA.y,'ro') % pax peaks in HF
    for k = 1:length(adhBound)
        adhBoundary = adhBound{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 1) %adhesion boundary
    end
    print('-depsc2', '-r300', strcat(epsPath,'/paxPeakDistribution',num2str(ii,iiformat),'.eps'));
    hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')
   
    close(h2)
    close(h1)
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