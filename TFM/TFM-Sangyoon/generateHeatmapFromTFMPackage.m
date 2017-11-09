function [strainEnergy,totalIntSelecChan,pixelTraction,pixelIntSelecChan] = generateHeatmapFromTFMPackage( pathForTheMovieDataFile,band,tmax,varargin )
%generateHeatmapFromTFMPackage generates heatmap from forcefield stored in
%movieData.
% input:    pathForTheMovieDataFile:    path to the movieData file
%           band:                       band width for cutting border
%           (default=4)
%           a certain point (default = false)
% output:   
%           strainEnergy: 1/2*integral(traction*u) in femtoJ (10^-15 Joule)
%           images of heatmap stored in pathForTheMovieDataFile/heatmap
ip =inputParser;
% ip.addRequired('pathForTheMovieDataFile',@ischar);%@(x)isscalar(x)||isempty(x));
% ip.addOptional('band',0,@isscalar);
% ip.addOptional('tmax',[],@isscalar);
ip.addParamValue('chanIntensity',[],@isnumeric); % channel to quantify intensity (2 or 3)
ip.addParamValue('vectorScale',1,@isnumeric); % channel to quantify intensity (2 or 3)
ip.addParamValue('tmin',[],@isnumeric); % channel to quantify intensity (2 or 3)
% ip.parse('pathForTheMovieDataFile','band','tmax',varargin{:});
ip.parse(varargin{:});
% pathForTheMovieDataFile=ip.Results.pathForTheMovieDataFile;
% band = ip.Results.band;
% tmax = ip.Results.tmax;
selectedChannel=ip.Results.chanIntensity;
vectorScale=ip.Results.vectorScale;
tmin=ip.Results.tmin;

if nargin < 2
    band = 4;
    tmax=[];
elseif nargin <3
    tmax=[];
end
% Load the MovieData
if isa(pathForTheMovieDataFile,'MovieData')
    movieData = pathForTheMovieDataFile;
%     movieData = MovieData.load(fullfile(MD.getPath,MD.movieDataFileName_));
    pathForTheMovieDataFile = movieData.getPath;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
    % structMD = load(movieDataPath);
    % movieData = structMD.MD;
    movieData = MovieData.load(movieDataPath);
end
% Get whole frame number
nFrames = movieData.nFrames_;
% Get TFM package
TFMPackage = movieData.getPackage(movieData.getPackageIndex('TFMPackage'));
% Load the displField
iDispFieldProc = 3;
displFieldProc=TFMPackage.processes_{iDispFieldProc};
if isempty(displFieldProc)
    iDispFieldProc = 2;
    displFieldProc=TFMPackage.processes_{iDispFieldProc};
end

% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
forceField=forceFieldProc.loadChannelOutput;

% filter forcefield temporally
% for each node
% if nFrames>1
%     Nnodes = length(forceField(1).pos);
%     for k=1:Nnodes
%         % see the profile
%         curVecX = arrayfun(@(x) x.vec(k,1),forceField);
%         curVecY = arrayfun(@(x) x.vec(k,2),forceField);
%     %     figure, plot(1:length(curVecX),curVecX)
%         [~,px] = csaps(1:length(curVecX),curVecX);
%         [~,py] = csaps(1:length(curVecY),curVecY);
%         smVecX2 = csaps(1:length(curVecX),curVecX,px*0.5);
%         smVecY2 = csaps(1:length(curVecY),curVecY,py*0.5);
%         smVecX3 = fnval(smVecX2,1:length(curVecX));
%         smVecY3 = fnval(smVecY2,1:length(curVecY));
%     %     hold all, plot(1:length(curVecX),smVecX3)
%         for ii=1:length(curVecX)
%             forceField(ii).vec(k,:) = [smVecX3(ii) smVecY3(ii)];
%         end
%     end
% end

% Set up the output file path
outputFilePath = [pathForTheMovieDataFile filesep 'Heatmaps'];
paxPath = [outputFilePath filesep 'pax'];
tifPath = [outputFilePath filesep 'tifs'];
figPath = [outputFilePath filesep 'figs'];
forcemapPath = [outputFilePath filesep 'fMap'];
epsPath = [outputFilePath filesep 'eps'];
if ~exist(tifPath,'dir') || ~exist(paxPath,'dir') || ~exist(epsPath,'dir') || ~exist(forcemapPath,'dir')
    mkdir(paxPath);
    mkdir(tifPath);
    mkdir(figPath);
    mkdir(forcemapPath);
    mkdir(epsPath);
end
% See if there is stage drift correction
iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=movieData.processes_{iSDCProc};
    if ~SDCProc.checkChannelOutput(1)
        error(['The channel must have been corrected ! ' ...
            'Please apply stage drift correction to all needed channels before '...
            'running displacement field calclation tracking!'])
    end
    if length(SDCProc.funParams_.ChannelIndex)>1
        iChan = 2;
    elseif length(SDCProc.funParams_.ChannelIndex) == 1
        iChan = SDCProc.funParams_.ChannelIndex;
    else
        error('No channel associated with SDC process!')
    end
    if iChan==2
        iBeadChan=1;
    else
        iBeadChan = SDCProc.funParams_.ChannelIndex(1);
    end
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
else
    iChan = 2;
end

iDisplFieldCalProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);
displFieldCalProc=movieData.processes_{iDisplFieldCalProc};
pDisp = parseProcessParams(displFieldCalProc);

% Use mask of first frame to filter displacementfield
if ~isempty(movieData.rois_) || ~isempty(movieData.roiMaskPath_)
%     maskArray = movieData.getROIMask;
    maskArray = imread(movieData.roiMaskPath_);
    if ~isempty(iSDCProc)

        %Parse input, store in parameter structure
        refFrame = double(imread(SDCProc.outFilePaths_{2,pDisp.ChannelIndex}));

        % Use mask of first frame to filter bead detection
        firstMask = refFrame>0; %false(size(refFrame));
        tempMask = maskArray(:,:,1);
        % firstMask(1:size(tempMask,1),1:size(tempMask,2)) = tempMask;
        tempMask2 = false(size(refFrame));
        y_shift = find(any(firstMask,2),1);
        x_shift = find(any(firstMask,1),1);
        tempMask2(y_shift:y_shift+size(tempMask,1)-1,x_shift:x_shift+size(tempMask,2)-1) = tempMask;
        firstMask = tempMask2 & firstMask;
        displFieldOriginal=displFieldProc.loadChannelOutput;
        displField = filterDisplacementField(displFieldOriginal,firstMask);
    else        
        firstMask = maskArray(:,:,1);
        displFieldOriginal=displFieldProc.loadChannelOutput;
        displField = filterDisplacementField(displFieldOriginal,firstMask);
    end
else
    displField=displFieldProc.loadChannelOutput;
end

%Find the maximum force.
[reg_grid1,~,~,~]=createRegGridFromDisplField(displField,1); %2=2 times fine interpolation

% band width for cutting border
%     band=4;
if isempty(tmax)
    display('Estimating maximum force magnitude ...')
    tic
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
    tmax = 0.8*tmax;
    display(['Estimated force maximum = ' num2str(tmax) ' Pa.'])
    toc
end
if isempty(tmin)
    tmin = 100000000;
    display('Estimating minimum force magnitude ...')
    tic
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
    toc
end    
%     tmax = 2590;
%     tmin = tmin-0.1;
%     tmax=tmax/5;
%     LeftUpperCorner(1:2) = [min(displField(1).pos(:,1)), min(displField(1).pos(:,2))];
%     RightLowerCorner(1:2) = [max(displField(1).pos(:,1)), max(displField(1).pos(:,2))];

[reg_grid,~,~,spacing]=createRegGridFromDisplField(displField,4); %2=2 times fine interpolation

hl = []; %handle for scale bar
iiformat = ['%.' '3' 'd'];
TSlevel = zeros(nFrames,1);
ii=1;
%     paxLevel = zeros(nFrames,1);

% Load Cell Segmentation
iMask = movieData.getProcessIndex('MaskRefinementProcess');
if ~isempty(iMask)
    maskProc = movieData.getProcess(iMask);
    bwPI4 = maskProc.loadChannelOutput(iChan,1);
    if ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        I = padarray(bwPI4, [maxY, maxX]);
        bwPI4 = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
    else
        iMask = movieData.getProcessIndex('ThresholdProcess');
        maskProc = movieData.getProcess(iMask);
        bwPI4 = maskProc.loadChannelOutput(iChan,ii);
        if ~isempty(iSDCProc)
            maxX = ceil(max(abs(T(:, 2))));
            maxY = ceil(max(abs(T(:, 1))));
            Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
            I = padarray(bwPI4, [maxY, maxX]);
            bwPI4 = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
        end
    end
else
    % if there was no cell mask, just use the entire pixel as a mask
    firstBeadImg=SDCProc.loadChannelOutput(iBeadChan,1);
    bwPI4 = true(size(firstBeadImg,1),size(firstBeadImg,2));
end
strainEnergy = zeros(nFrames,1);
% Boundary cutting - I'll take care of this boundary effect later
if band>0
    reg_grid(1:band,:,:)=[];
    reg_grid(:,1:band,:)=[];
    reg_grid(end-band+1:end,:,:)=[];
    reg_grid(:,end-band+1:end,:)=[];
end

maskCrop = bwPI4(reg_grid(1,1,2):reg_grid(end,end,2),reg_grid(1,1,1):reg_grid(end,end,1));
nSegPixel = sum(maskCrop(:));

pixelTraction = zeros(nSegPixel,nFrames);
if ~isempty(selectedChannel)
    totalIntSelecChan = zeros(nFrames,1);
    pixelIntSelecChan = zeros(nSegPixel,nFrames);
end
h1 = figure('color','w');
for ii=1:nFrames
    [grid_mat,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
    pos = [reshape(grid_mat(:,:,1),[],1) reshape(grid_mat(:,:,2),[],1)]; %dense
    disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    [~,if_mat,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid);
    force_vec = [reshape(if_mat(:,:,1),[],1) reshape(if_mat(:,:,2),[],1)]; 

    [~,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],grid_mat); %1:cluster size
    tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;

    grid_mat_quiver=grid_mat(3:end-2,3:end-2,:);
    
    % drawing
%         hs = surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm,'FaceColor','interp',...
%             'EdgeColor','none', 'FaceLighting','gouraud');%, 'FaceLighting','phong');
%         zlim([tmin tmax]), view(0,90)
%     hs = pcolor(grid_mat(:,:,1), grid_mat(:,:,2), tnorm);%,[tmin tmax]);
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1);
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2);
    set(h1, 'Position', [100 100 (imSizeX+1)*1.25 imSizeY+1])
    
    subplot(1,1,1)
    subplot('Position',[0 0 0.8 1])
%     hs = pcolor(grid_mat(:,:,1), grid_mat(:,:,2), tnorm);%,[tmin tmax]);
%     colormap jet;
%     shading interp
%     caxis([tmin tmax])
%     set(gca, 'DataAspectRatio', [1,1,1],'Ydir','reverse');
    [XI,YI]=meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
    tsMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tnorm,XI,YI,'cubic');
    tsMapX = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tmat(:,:,1),XI,YI,'cubic');
    tsMapY = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tmat(:,:,2),XI,YI,'cubic');
    uMapX = griddata(grid_mat(:,:,1),grid_mat(:,:,2),iu_mat(:,:,1),XI,YI,'cubic');
    uMapY = griddata(grid_mat(:,:,1),grid_mat(:,:,2),iu_mat(:,:,1),XI,YI,'cubic');
    imshow(tsMap,[tmin tmax]), colormap jet;
    % strain energy calculation
    curStrainEnergy = 1/2*maskCrop.*(tsMapX.*uMapX+tsMapY.*uMapY)...
                                        *(movieData.pixelSize_*1e-9)^3/1e-15; %femto Joule
    strainEnergy(ii)= sum(curStrainEnergy(:)); % fJ
    pixelTraction = tsMap(maskCrop(:));

    % unit vector plot
    hold on
    [reg_grid_coarse,~,~,spacing]=createRegGridFromDisplField(displField,0.5); %2=2 times fine interpolation
    [grid_mat_coarse,iu_mat_coarse,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid_coarse);
    pos_coarse = [reshape(grid_mat_coarse(:,:,1),[],1) reshape(grid_mat_coarse(:,:,2),[],1)]; %dense
    disp_vec_coarse = [reshape(iu_mat_coarse(:,:,1),[],1) reshape(iu_mat_coarse(:,:,2),[],1)]; 
    [~,if_mat_coarse,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid_coarse);
    force_vec_coarse = [reshape(if_mat_coarse(:,:,1),[],1) reshape(if_mat_coarse(:,:,2),[],1)]; 

    [~,tmat_coarse, ~, ~] = interp_vec2grid(pos_coarse+disp_vec_coarse, force_vec_coarse,[],grid_mat_coarse); %1:cluster size

%     tmat_coarse(end-band/4-1:end,:,:)=[];
%     tmat_coarse(:,end-band/4-1:end,:)=[];
%     tmat_coarse(1:1+band/4+1,:,:)=[];
%     tmat_coarse(:,1:1+band/4+1,:)=[];
%     grid_mat_coarse(end-band/4-1:end,:,:)=[];
%     grid_mat_coarse(:,end-band/4-1:end,:)=[];
%     grid_mat_coarse(1:1+band/4+1,:,:)=[];
%     grid_mat_coarse(:,1:1+band/4+1,:)=[];

    tmat_vecx = reshape(tmat_coarse(:,:,1),[],1);
    tmat_vecy = reshape(tmat_coarse(:,:,2),[],1);
    pos_vecx = reshape(grid_mat_coarse(:,:,1),[],1);
    pos_vecy = reshape(grid_mat_coarse(:,:,2),[],1);
    forceScale=0.1*max(sqrt(tmat_vecx.^2+tmat_vecy.^2));
%     hq = quiver(pos_vecx,pos_vecy, tmat_vecx./forceScale,tmat_vecy./forceScale,0,'k');
    hq = quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), vectorScale*tmat_vecx./forceScale,vectorScale*tmat_vecy./forceScale,0,'Color',[75/255 0/255 130/255]);
    hq.ShowArrowHead = 'off';
%     hq.LineWidth=0.5;
%     anno=hq.Annotation;
%     anno.LegendInformation
%     hq.AutoScale='off';
%     hq = arrayfun(@(x,y,u,v) annotation('arrow',[(x-grid_mat(1,1,1))/(imSizeX+1),(x-grid_mat(1,1,1)+vectorScale*u./forceScale)/(imSizeX+1)],...
%         [(y-grid_mat(1,1,2))/(imSizeY+1),(y-grid_mat(1,1,2)+vectorScale*v./forceScale)/(imSizeY+1)],'Units','points'),...
%         pos_vecx,pos_vecy,tmat_vecx,tmat_vecy);
%     hq.HeadStyle = 'cback2';
%     hq.HeadLength = 15; % does not work!
    % boundary
    bwPI4 = maskProc.loadChannelOutput(iChan,ii);
    if ~isempty(iSDCProc)
        maxX = ceil(max(abs(T(:, 2))));
        maxY = ceil(max(abs(T(:, 1))));
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
        I = padarray(bwPI4, [maxY, maxX]);
        bwPI4 = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
    end
    
    maskCrop = bwPI4(reg_grid(1,1,2):reg_grid(end,end,2),reg_grid(1,1,1):reg_grid(end,end,1));
    [cB,~,nCBD]  = bwboundaries(maskCrop,'noholes');
    for kk=1:nCBD
        boundary = cB{kk};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
    end

    % Scale bar 2000nm
%     if isempty(hl)
%         hold on
    scale = 5; %micron
    hl = line([10 10+round(scale*1000/movieData.pixelSize_)],[15 15],'LineWidth',2,'Color',[1,1,1]);
    text(10,30,[num2str(scale), ' um.'],'Color','w');
%     disp(['Scale bar: ', num2str(scale), ' um.'])
%     end
    axis off
    hold off
    subplot('Position',[0.8 0.1 0.1 0.8])
    axis tight
    caxis([tmin tmax]), axis off
%     if isempty(hc)
    hc = colorbar('West');
    set(hc,'Fontsize',18)
%     end
    nChannels  = length(movieData.channels_);
    if nChannels==2
        % loading paxillin image
        if ~isempty(iSDCProc)
            secondImage=(SDCProc.loadChannelOutput(iChan,ii)); %movieData.channels_(2).loadImage(ii);
        else
            secondImage=movieData.getChannel(iChan).loadImage(ii); 
        end
    %     paxImageCropped = paxImage(indULy+spacing*band:indBRy-spacing*band,indULx+spacing*band:indBRx-spacing*band);
        secondImageCropped = secondImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
        %Scale bar
        secondImageCropped(15:16,10:10+round(2000/movieData.pixelSize_))=max(max(secondImageCropped)); % this is 2 um
        imwrite(secondImageCropped, strcat(paxPath,'/paxCroppedTif',num2str(ii,iiformat),'.tif'));
        % composite for both channels
        compImage(:,:,1) = imadjust(tsMap/tmax,[],[]);
        doubleSecondImg = double(secondImageCropped)/double(max(secondImageCropped(:)));
        if ~isempty(iSDCProc)
            secondImageUnshifted=movieData.getChannel(iChan).loadImage(ii); 
            doubleSecondImgUnshifted = double(secondImageUnshifted)/double(max(secondImageCropped(:)));
            minPax= min(doubleSecondImgUnshifted(:));
        else
            minPax= min(doubleSecondImg(:));
        end
        
        compImage(:,:,2) = imadjust(doubleSecondImg,[minPax,max(doubleSecondImg(:))],[]);
        compImage(:,:,3) = imadjust(tsMap/tmax,[],[]);
        imwrite(compImage, strcat(paxPath,'/CombPaxForceTif',num2str(ii,iiformat),'.tif'));
%         figure, imshow(compImage,[])
    elseif nChannels==3
        % loading paxillin image
        if ~isempty(iSDCProc)
            secondImage=(SDCProc.loadChannelOutput(iChan,ii)); %movieData.channels_(2).loadImage(ii);
            thirdImage=(SDCProc.loadChannelOutput(iChan+1,ii)); %movieData.channels_(2).loadImage(ii);
        else
            secondImage=movieData.getChannel(iChan).loadImage(ii); 
            thirdImage=movieData.getChannel(iChan+1).loadImage(ii); 
        end
        secondImageCropped = secondImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
        thirdImageCropped = thirdImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
        %Scale bar
        secondImageCropped(15:16,10:10+round(2000/movieData.pixelSize_))=max(max(secondImageCropped)); % this is 2 um
        thirdPath = [outputFilePath filesep 'thirdChannel'];
        if ~exist(thirdPath,'dir') 
            mkdir(thirdPath);
        end
        thirdImageCropped(15:16,10:10+round(2000/movieData.pixelSize_))=max(max(thirdImageCropped)); % this is 2 um
        imwrite(secondImageCropped, strcat(paxPath,'/secondCroppedTif',num2str(ii,iiformat),'.tif'));
        imwrite(thirdImageCropped, strcat(thirdPath,'/thirdCroppedTif',num2str(ii,iiformat),'.tif'));
        
        h2 = figure('color','w');
        set(h2, 'Position', [100 100 (imSizeX+1) imSizeY+1])

        thirdImageUnshifted=movieData.getChannel(iChan+1).loadImage(ii); 
        doubleThirdImageUnshifted = double(thirdImageUnshifted)/double(max(thirdImageUnshifted(:)));
        thirdImageLog = log(double(doubleThirdImageUnshifted));
        thirdImageInverted=(ones(size(thirdImageLog)))*max(thirdImageLog(:))-thirdImageLog;
        minThird= max(min(thirdImageInverted(:)),0);
        maxThird= max(thirdImageInverted(:));
        
        thirdImageLog = log(double(thirdImageCropped));
        thirdImageInverted=(ones(size(thirdImageLog)))*max(thirdImageLog(:))-thirdImageLog;
        
%         imshow(thirdImageInverted,[minThird+0.3*maxThird maxThird]), hold on
        imshow(thirdImageInverted,[minThird maxThird])
        hold on
%         colormap jet;
        for kk=1:nCBD
            boundary = cB{kk};
            plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1) % cell boundary
        end
        % saving
        I2 = getframe(h2);
        imwrite(I2.cdata, strcat(tifPath,'/thirdImgTif',num2str(ii,iiformat),'.tif'));
        hgsave(h2,strcat(figPath,'/thirdImgFig',num2str(ii,iiformat)),'-v7.3')
        print(h2,strcat(epsPath,'/thirdImgEps',num2str(ii,iiformat),'.eps'),'-depsc2')
        hold off
        
        % average intensity quantification on selected image
        if ~isempty(selectedChannel)
            if selectedChannel==2
                curIntSelecChan = 1/2*double(bwPI4).*double(secondImage); %AU
                totalIntSelecChan(ii)= sum(curIntSelecChan(:)); % AU
                pixelIntSelecChan(:,ii) = secondImage(maskCrop(:));
                pixelID = find(maskCrop); %pixel id
            elseif selectedChannel==3
                curIntSelecChan = 1/2*double(bwPI4).*double(thirdImage); %AU
                totalIntSelecChan(ii)= sum(curIntSelecChan(:)); % fJ
                pixelIntSelecChan(:,ii) = thirdImageCropped(maskCrop(:));
                pixelID = find(maskCrop); %pixel id
                % Showing plot between pixelTraction and pixelIntSelecChan
                hScatter = figure; plot(pixelIntSelecChan(:,ii),pixelTraction(:,ii),'.')
                % Ask limit for high traction and high vim
                highTraction = input('Limit for high traction above which you want to plot? :');
                highVim = input('Limit for high vimentin level above which you want to plot? :');
                % Showing these regions by boundaries
                indHighTraction = pixelTraction(:,ii)>highTraction;
                indHighVim = pixelIntSelecChan(:,ii)>highVim;
                hold on
                plot(pixelIntSelecChan(indHighTraction,ii),pixelTraction(indHighTraction,ii),'r.')
                plot(pixelIntSelecChan(indHighVim,ii),pixelTraction(indHighVim,ii),'g.')
                close(hScatter);
                %Showing them in 2D histogram)
                if max(pixelTraction(:,ii))<100
                    yBins = round(min(pixelTraction(:,ii))):round(max(pixelTraction(:,ii)));
                    xBins = min(pixelIntSelecChan(:,ii)):1:max(pixelIntSelecChan(:,ii));
                else
                    yBins = linspace(round(min(pixelTraction(:,ii))),round(max(pixelTraction(:,ii))),100);
                    xBins = linspace(min(pixelIntSelecChan(:,ii)),max(pixelIntSelecChan(:,ii)),100);
                end
    %             yBins = round(min(pixelTraction)):100:round(max(pixelTraction));
                hHist2D = figure; hold on
                densityplot(pixelIntSelecChan(:,ii), pixelTraction(:,ii), xBins, yBins,'DisplayFunction', @log);
                h_cb=colorbar;
                h_cb.Label.String = 'Occurence, 10 ^';
                ax = gca;
                axpos = ax.Position;
                cpos = h_cb.Position;
                cpos(3) = 0.5*cpos(3);
                cpos(2) = cpos(2)+0.05*cpos(4);
                cpos(4) = 0.9*cpos(4);
                h_cb.Position = cpos;
                ax.Position = axpos;

                ylabel('Traction (Pa)')
                xlabel('Vimentin Intensity (A.U.)')
                % rectacgle
                if sum(indHighTraction)>5
                    rectangle('Position',[min(pixelIntSelecChan(indHighTraction,ii)) highTraction ...
                        max(pixelIntSelecChan(indHighTraction,ii))-min(pixelIntSelecChan(indHighTraction,ii)) ...
                        max(pixelTraction(indHighTraction,ii))-highTraction],'EdgeColor','r')
                end
                if sum(indHighVim)>5
                    rectangle('Position',[highVim min(pixelTraction(indHighVim,ii)) ...
                        max(pixelIntSelecChan(indHighVim,ii))-highVim ...
                        max(pixelTraction(indHighVim,ii))-min(pixelTraction(indHighVim,ii))],'EdgeColor','g')
                end
                % save
                print('-depsc2', '-r150', strcat(epsPath,'/Hist2DbtwVimAndTraction',num2str(ii,iiformat),'.eps'));
                close(hHist2D)

    %             map = getScatterQuantification(pixelIntSelecChan,pixelTraction,xBins,yBins);
    %             figure, imshow(map,[0 0.0002]), colormap jet

                % by making mask
                highTracMask = false(size(maskCrop));
                highTracMask(pixelID(indHighTraction)) = true;
                highVimMask = false(size(maskCrop));
                highVimMask(pixelID(indHighVim)) = true;
                % and by making it boundaries
                [tB,~,nTBD]  = bwboundaries(highTracMask,'noholes');
                hVim=figure; imshow(thirdImageCropped,[]), hold on
                for kk=1:nTBD
                    boundary = tB{kk};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5) % high traction boundary
                end
                % cell mask
                [cB,~,nCBD]  = bwboundaries(maskCrop,'noholes');
                for kk=1:nCBD
                    boundary = cB{kk};
                    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 0.5) % cell boundary
                end
                print('-depsc2', '-r150', strcat(epsPath,'/thridImageWithHighTraction',num2str(ii,iiformat),'.eps'));
                close(hVim)

                [vB,~,nVBD]  = bwboundaries(highVimMask,'noholes');
                hT=figure; imshow(tsMap,[tmin tmax]), colormap jet, hold on
                for kk=1:nVBD
                    boundary = vB{kk};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5) % cell boundary
                end
                for kk=1:nCBD
                    boundary = cB{kk};
                    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 0.5) % cell boundary
                end
                print('-depsc2', '-r150', strcat(epsPath,'/tractionImageWithHighVim',num2str(ii,iiformat),'.eps'));
                close(hT)

            end       
            % composite for both channels
            compImage(:,:,1) = imadjust(tsMap/tmax,[],[]);
            doubleSecondImg = double(secondImageCropped)/double(max(secondImageCropped(:)));
            doubleThirdImg = double(thirdImageCropped)/double(max(thirdImageCropped(:)));
            if ~isempty(iSDCProc)
                secondImageUnshifted=movieData.getChannel(iChan).loadImage(ii); 
                doubleSecondImgUnshifted = double(secondImageUnshifted)/double(max(secondImageCropped(:)));
                minPax= min(doubleSecondImgUnshifted(:));
                thirdImageUnshifted=movieData.getChannel(iChan+1).loadImage(ii); 
                doubleThirdImageUnshifted = double(thirdImageUnshifted)/double(max(thirdImageUnshifted(:)));
                minThird= min(doubleThirdImageUnshifted(:));
            else
                minPax= min(doubleSecondImg(:));
                minThird= min(doubleThirdImg(:));
            end

            compImage(:,:,3) = imadjust(doubleSecondImg,[minPax,max(doubleSecondImg(:))],[]);
            compImage(:,:,2) = imadjust(doubleThirdImg,[minThird,max(doubleThirdImg(:))],[]);
            imwrite(compImage, strcat(paxPath,'/CombPaxForceTif',num2str(ii,iiformat),'.tif'));
    %         figure, imshow(compImage,[])
        end
    end

    % saving
    I = getframe(h1);
    imwrite(I.cdata, strcat(tifPath,'/stressMagTif',num2str(ii,iiformat),'.tif'));
    imwrite(uint16(round(tsMap*10)),strcat(forcemapPath,'/force',num2str(ii,iiformat),' divide by 10 for correct mag','.tif'));

%         hgexport(h1,strcat(tifPath,'/stressMagTif',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h1,strcat(figPath,'/stressMagFig',num2str(ii,iiformat)),'-v7.3')

%     set(h1,'Renderer','zbuffer')
%     saveas(h1,strcat(epsPath,'/stressMagEps',num2str(ii,iiformat),'.eps'),'eps') 
    print(h1,strcat(epsPath,'/stressMagEps',num2str(ii,iiformat),'.eps'),'-depsc2')
    hold off
%     delete(hs)
%     delete(hq)
%     delete(hl);
%     delete(hc);
%     hl = []; %handle for scale bar
%     
end
close(h1)
if nChannels>2
    close(h2)
end
return;
% to run the function:
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret/TFM/cell 5/c647_im',6);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret/TFM/cell3/TFM',40);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130110 cell 4 cropped',4);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130307 data/1301331 Cell3_pax TIRF',8);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF',4);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF only NA',4);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130410 paxilin crop/Cell12/crop_movieData',8,false,1600);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130110 Cell7 2Frames/ROIAnalysis',8,false, 2100);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130410 paxilin crop/Cell12/crop_movieData/ROIAnalysis',16,false,2500)
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Youbean/130131 cell3paxTIRF 2Frames/ROIAnalysis',16,false,4000)
% desktop version
generateHeatmapFromTFMPackage('/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret/TFM/cell 5/c647_im',6);
generateHeatmapFromTFMPackage('/home/sh268/files/LCCB/shared/X-change/forSangyoon/fromYoubean/120907 Cell5/',4);