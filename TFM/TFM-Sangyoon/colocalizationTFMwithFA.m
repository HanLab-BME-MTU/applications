function [] = colocalizationTFMwithFA( pathForTheMovieDataFile,outputPath,band,pointTF,tmax )
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
elseif nargin <4
    pointTF = false;
    tmax=[];
elseif nargin<5
    tmax=[];
end

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
%     tmax = 2590;
%     tmin = tmin-0.1;
%     tmax=tmax/5;
%     LeftUpperCorner(1:2) = [min(displField(1).pos(:,1)), min(displField(1).pos(:,2))];
%     RightLowerCorner(1:2) = [max(displField(1).pos(:,1)), max(displField(1).pos(:,2))];

[reg_grid,~,~,spacing]=createRegGridFromDisplField(displField,4); %4=4 times fine interpolation

%     h2 = figure;
hold off
%     set(h2, 'Position', [100+imSizeX*10/9 100 imSizeX imSizeY])
hc = []; %handle for colorbar
hl = []; %handle for scale bar
iiformat = ['%.' '3' 'd'];
TSlevel = zeros(nFrames,1);
%     paxLevel = zeros(nFrames,1);
SegmentationPackage = movieData.getPackage(movieData.getPackageIndex('SegmentationPackage'));
minSize = round((1000/movieData.pixelSize_)*(500/movieData.pixelSize_)); %adhesion limit=1um*.5um

forcePeaksInfo(nFrames,1)=struct('allfPeaks',[],'indInFAN',[],...
    'indInNAN',[],'nAll',[],'nInFAN',[],'nInNAN',[]);
paxPeaksInfo(nFrames,1)=struct('allpPeaks',[],'indInHFN',[],...
    'indInLFN',[],'nAll',[],'nInHFN',[],'nInLFN',[]);

for ii=1:nFrames
    [grid_mat,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
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
    [~,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],grid_mat); %1:cluster size
    tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;

    % output to a file
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1);
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2);

    [XI,YI]=meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
    tsMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tnorm,XI,YI,'cubic');
    % Mask 
    maskProc = SegmentationPackage.processes_{2};
    mask = maskProc.loadChannelOutput(2,ii);
    cropMask = mask(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
    % Boundary
    B  = bwboundaries(cropMask,'noholes');
    boundary = B{1};

    % Obtain peaks of traction stress
    sigmaPSF = 2;
    alpha = 0.01;
    pstruct = pointSourceDetection(tsMap, sigmaPSF, 'alpha', alpha,'Mask',cropMask);
    
    % Pick up points in background:
    backgroundMask = imcomplement(cropMask);
    pstruct_BG = pointSourceDetection(tsMap, sigmaPSF, 'alpha', alpha,'Mask',backgroundMask);
    
    % Find 90 percentile force level of the peaks in background 
    meanA_BG = prctile(pstruct_BG.A,90);

    % Filter out peaks in the cell that are lower than meanA_BG
    ind = pstruct.A>meanA_BG;
    pstruct_filtered.x = pstruct.x(ind);
    pstruct_filtered.y = pstruct.y(ind);
    pstruct_filtered.A = pstruct.A(ind);
    
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
    maskAdhesion = blobSegmentThreshold(paxImageCropped,minSize,0,cropMask);
    adhBound = bwboundaries(maskAdhesion,'noholes');
    
    % force peaks in neighborhood of FAs
    indInFA = maskVectors(pstruct_filtered.x,pstruct_filtered.y,maskAdhesion);
    % with some neighborhood (5 pix)
    maskAdhesionN = bwdist(maskAdhesion)<=5;
    indInFAN = maskVectors(pstruct_filtered.x,pstruct_filtered.y,maskAdhesionN);
    
    % Nascent adhesion detection
    sigmaPSF_NA = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, 601*1e-9);
    naMask = ~maskAdhesion & cropMask;
    pstruct_NA = pointSourceDetection(paxImageCropped, sigmaPSF_NA, 'alpha', 0.05,'Mask',naMask);
%     plot(pstruct_NA.x,pstruct_NA.y,'co')

    % Nascent adhesion mask assuming isotropic shape (circles)
    maskNA = false(size(maskAdhesion));
    Nna = length(pstruct_NA.x);
    for k=1:Nna
        maskNA(round(pstruct_NA.y(k)),round(pstruct_NA.x(k))) = true;
    end
    maskNAN = bwdist(maskNA)<=500/movieData.pixelSize_; %500nm
    % Force peaks in neighborhood of FAs
    indInNAN = maskVectors(pstruct_filtered.x,pstruct_filtered.y,maskNAN);
    
    h1 = figure('color','w');
    set(h1, 'Position', [100 100 imSizeX*1.25 imSizeY])
    subplot('Position',[0 0 0.8 1])
    imshow(tsMap,[tmin tmax]), colormap gray;
    hold on;
%     plot(pstruct.x,pstruct.y,'mo') % all peaks
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
%     for k = 1:length(adhBound)
%         adhBoundary = adhBound{k};
%         plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 1) %adhesion boundary
%     end
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
    % Scale bar 2000nm
    plot(pstruct_filtered.x(indInFAN),pstruct_filtered.y(indInFAN),'ro') % force peak in FA
    plot(pstruct_filtered.x(indInNAN),pstruct_filtered.y(indInNAN),'yo') % force peak in NA
    plot(pstruct_filtered.x(~indInFAN&~indInNAN),pstruct_filtered.y(~indInFAN&~indInNAN),'go') % all peaks
    
    if isempty(hl)
        hold on
        hl = line([grid_mat(2,2,1),grid_mat(2,2,1)+2000/movieData.pixelSize_],[grid_mat(2,2,2),grid_mat(2,2,2)],'Color','w','Linewidth',2);
    end
    
    subplot('Position',[0.8 0.1 0.1 0.8])
    axis tight
    caxis([tmin tmax]), axis off
    hc = colorbar('West');
    set(hc,'Fontsize',16)

    print('-depsc2', '-r300', strcat(epsPath,'/forcePeakDistribution',num2str(ii,iiformat),'.eps'));
    hgexport(h1,strcat(forcetifPath,'/forcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h1,strcat(figPath,'/forcePeakFig',num2str(ii,iiformat)),'-v7.3')
    
    % Plot the number of force peaks 
%     fPeaks(1) = sum(numel(pstruct_filtered.x));
%     fPeaks(2) = sum(indInFAN);
%     fPeaks(3) = sum(indInNAN);
%     figure,bar(fPeaks)
%     print('-depsc2', '-r300', strcat(epsPath,'/forcePeakBar',num2str(ii,iiformat),'.eps'));
    forcePeaksInfo(ii).allfPeaks = pstruct_filtered;
    forcePeaksInfo(ii).indInFAN = indInFAN;
    forcePeaksInfo(ii).indInNAN = indInNAN;
    forcePeaksInfo(ii).nAll = length(pstruct_filtered.x);
    forcePeaksInfo(ii).nInFAN = sum(indInFAN);
    forcePeaksInfo(ii).nInNAN = sum(indInNAN);
    % now it's time to classify paxPeaks
    % Get all peaks
    pstruct_FA = pointSourceDetection(paxImageCropped, sigmaPSF_NA*1.6, 'alpha', alpha,'Mask',cropMask);
    % Mask for high force (bigger than background force)
    highForce = tmax*2/5;
    maskHighForce = im2bw(tsMap/max(max(tsMap)),highForce/max(max(tsMap)));
    indPaxInHF = maskVectors(pstruct_FA.x,pstruct_FA.y,maskHighForce);
    % Pax peaks in point force neighborhood whose level is lower than high
    % force leve
    % Filter out peaks in the cell that are lower than meanA_BG
    indLowPeaks = pstruct.A<highForce;
    pstruct_lowPeaks.x = pstruct.x(indLowPeaks);
    pstruct_lowPeaks.y = pstruct.y(indLowPeaks);
    pstruct_lowPeaks.A = pstruct.A(indLowPeaks);
    % low force mask assuming isotropic shape (circles)
    maskLowForce = false(size(paxImageCropped));
    N_LF = length(pstruct_lowPeaks.x);
    for k=1:N_LF
        maskLowForce(round(pstruct_lowPeaks.y(k)),round(pstruct_lowPeaks.x(k))) = true;
    end
    maskLowForceN = bwdist(maskLowForce)<=10;
    % Force peaks in neighborhood of FAs
    indPaxInLFN = maskVectors(pstruct_FA.x,pstruct_FA.y,maskLowForceN);
    
    %Scale bar 2 um
    paxImageCropped(15:16,10:10+round(2000/movieData.pixelSize_))=max(max(paxImageCropped));
    h2=figure; imshow(paxImageCropped,[])
    hold on;
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
    plot(pstruct_FA.x(indPaxInHF),pstruct_FA.y(indPaxInHF),'ro') % pax peaks in HF
    plot(pstruct_FA.x(indPaxInLFN),pstruct_FA.y(indPaxInLFN),'yo') % pax peaks near low force peaks
    plot(pstruct_FA.x(~indPaxInHF & ~indPaxInLFN),pstruct_FA.y(~indPaxInHF & ~indPaxInLFN),'go') % pax peaks in rest
    for k = 1:length(adhBound)
        adhBoundary = adhBound{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'w', 'LineWidth', 1)
    end
    print('-depsc2', '-r300', strcat(epsPath,'/paxPeakDistribution',num2str(ii,iiformat),'.eps'));
    hgexport(h2,strcat(paxtifPath,'/paxWithForcePeak',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h2,strcat(figPath,'/paxPeakFig',num2str(ii,iiformat)),'-v7.3')

%     % Plot the number of paxillin peaks 
%     paxPeaks(1) = sum(numel(pstruct_FA.x));
%     paxPeaks(2) = sum(indPaxInHF);
%     figure,bar(paxPeaks)
%     print('-depsc2', '-r300', strcat(epsPath,'/paxPeakBar',num2str(ii,iiformat),'.eps'));
    paxPeaksInfo(ii).allpPeaks = pstruct_FA;
    paxPeaksInfo(ii).indInHFN = indPaxInHF;
    paxPeaksInfo(ii).indInLFN = indPaxInLFN;
    paxPeaksInfo(ii).nAll = length(pstruct_FA.x);
    paxPeaksInfo(ii).nInHFN = sum(indPaxInHF);
    paxPeaksInfo(ii).nInLFN = sum(indPaxInLFN);

%     % Get intensities at pstruct
%     Npoints = length(pstruct_filtered.x);
%     intenPax = zeros(Npoints,1);
%     magForce = zeros(Npoints,1);
%     for k=1:Npoints
%         intenPax(k) = paxImageCropped(round(pstruct_filtered.y(k)),round(pstruct_filtered.x(k)));
%         magForce(k) = tsMap(round(pstruct_filtered.y(k)),round(pstruct_filtered.x(k)));
%     end
%     h3=figure; plot(magForce,intenPax,'k.')
%     hgexport(h3,strcat(tifPath,'/forceVsPaxInten',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
   
    imwrite(tsMap, strcat(forcemapPath,'/TSTif',num2str(ii,iiformat),'.tif'));
    imwrite(paxImageCropped, strcat(paxPath,'/paxCroppedTif',num2str(ii,iiformat),'.tif'));
    close(h2)
    close(h1)
end
if pointTF
    t = (0:nFrames-1)*movieData.timeInterval_;
    figure,[AX,H1,H2] = plotyy(t,TSlevel,t,paxLevel,'plot');
    set(get(AX(1),'Ylabel'),'String','Traction Stress (Pa)') 
    set(get(AX(2),'Ylabel'),'String','Paxillin Fluorescence Intensity (A.U)') 
    xlabel('Time (sec)') 
    set(H1,'LineStyle','--')
end
save([dataPath filesep 'PeaksInfo.mat'],'forcePeaksInfo','paxPeaksInfo');
return;
% to run the function:
colocalizationTFMwithFA('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130110 Cell7 2Frames/ROIAnalysis',8,false);
colocalizationTFMwithFA('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis','L1 2nd',16,false,4000);
colocalizationTFMwithFA('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis','L2 0th',16,false,4000);