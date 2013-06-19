function [] = generateHeatmapFromTFMPackage( pathForTheMovieDataFile,band,pointTF,tmax )
%generateHeatmapFromTFMPackage generates heatmap from forcefield stored in
%movieData.
% input:    pathForTheMovieDataFile:    path to the movieData file
%           band:                       band width for cutting border
%           (default=4)
%           pointTF:                    true if you want to trace force at
%           a certain point (default = false)
% output:   images of heatmap stored in pathForTheMovieDataFile/heatmap
if nargin < 2
    band = 4;
    pointTF = false;
    tmax=[];
elseif nargin <3
    pointTF = false;
    tmax=[];
elseif nargin<4
    tmax=[];
end
% Load the MovieData
movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
% movieData = MovieData.load(movieDataPath);
load(movieDataPath);
movieData = MD;
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
outputFilePath = [pathForTheMovieDataFile filesep 'Heatmaps'];
paxPath = [pathForTheMovieDataFile filesep 'pax'];
tifPath = [outputFilePath filesep 'tifs'];
figPath = [outputFilePath filesep 'figs'];
epsPath = [outputFilePath filesep 'eps'];
if ~exist(tifPath,'dir') || ~exist(paxPath,'dir') || ~exist(epsPath,'dir')
    mkdir(paxPath);
    mkdir(tifPath);
    mkdir(figPath);
    mkdir(epsPath);
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

[reg_grid,~,~,spacing]=createRegGridFromDisplField(displField,4); %2=2 times fine interpolation

h1 = figure('color','w');
%     h2 = figure;
hold off
% [indULy,indULx] = ind2sub(size(firstMask),find(firstMask,1,'first'));
% [indBRy,indBRx] = ind2sub(size(firstMask),find(firstMask,1,'last'));
% imSizeX = indBRx-indULx-band*spacing;
% imSizeY = indBRy-indULy-band*spacing;

%     set(h2, 'Position', [100+imSizeX*10/9 100 imSizeX imSizeY])
hc = []; %handle for colorbar
hl = []; %handle for scale bar
iiformat = ['%.' '3' 'd'];
TSlevel = zeros(nFrames,1);
%     paxLevel = zeros(nFrames,1);

for ii=1:nFrames
    [grid_mat,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
    pos = [reshape(grid_mat(:,:,1),[],1) reshape(grid_mat(:,:,2),[],1)]; %dense
    disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    [~,if_mat,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid);
    force_vec = [reshape(if_mat(:,:,1),[],1) reshape(if_mat(:,:,2),[],1)]; 

    [~,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],grid_mat); %1:cluster size
    tnorm = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;

    % Boundary cutting - I'll take care of this boundary effect later
    tnorm(end-band:end,:)=[];
    tnorm(:,end-band:end)=[];
    tnorm(1:1+band,:)=[];
    tnorm(:,1:1+band)=[];
%     tmat(end-band-2:end,:,:)=[];
%     tmat(:,end-band-2:end,:)=[];
%     tmat(1:1+band+2,:,:)=[];
%     tmat(:,1:1+band+2,:)=[];
    grid_mat(end-band:end,:,:)=[];
    grid_mat(:,end-band:end,:)=[];
    grid_mat(1:1+band,:,:)=[];
    grid_mat(:,1:1+band,:)=[];

    grid_mat_quiver=grid_mat(3:end-2,3:end-2,:);
    
    % drawing
%         hs = surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm,'FaceColor','interp',...
%             'EdgeColor','none', 'FaceLighting','gouraud');%, 'FaceLighting','phong');
%         zlim([tmin tmax]), view(0,90)
%     hs = pcolor(grid_mat(:,:,1), grid_mat(:,:,2), tnorm);%,[tmin tmax]);
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1);
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2);
    if ii==1
        set(h1, 'Position', [100 100 imSizeX*1.25 imSizeY])
    end
    subplot(1,1,1)
    subplot('Position',[0 0 0.8 1])
%     hs = pcolor(grid_mat(:,:,1), grid_mat(:,:,2), tnorm);%,[tmin tmax]);
%     colormap jet;
%     shading interp
%     caxis([tmin tmax])
%     set(gca, 'DataAspectRatio', [1,1,1],'Ydir','reverse');
    [XI,YI]=meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
    tsMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),tnorm,XI,YI,'cubic');
    imshow(tsMap,[tmin tmax]), colormap jet;

    % unit vector plot
    hold on
    [reg_grid_coarse,~,~,spacing]=createRegGridFromDisplField(displField,1); %2=2 times fine interpolation
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
    forceScale=0.25*max(sqrt(tmat_vecx.^2+tmat_vecy.^2));
%     hq = quiver(pos_vecx,pos_vecy, tmat_vecx./forceScale,tmat_vecy./forceScale,0,'k');
    hq = quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), tmat_vecx./forceScale,tmat_vecy./forceScale,0,'Color',[75/255 0/255 130/255]);

    % Scale bar 2000nm
    if isempty(hl)
        hold on
        hl = line([grid_mat(2,2,1),grid_mat(2,2,1)+2000/movieData.pixelSize_],[grid_mat(2,2,2),grid_mat(2,2,2)],'Color','w','Linewidth',2);
    end
%         xlim([LeftUpperCorner(1) RightLowerCorner(1)])
%         ylim([LeftUpperCorner(2) RightLowerCorner(2)])
    % pick a point
    if pointTF && ii==1
        point = ginput(1);
    end
    % point = [357.7008  319.1465]
    % grid_mat's sub for point
%     spacing = grid_mat(2,1,1)-grid_mat(1,1,1);
    if pointTF
        indTS = find(grid_mat(:,:,1)>point(1)-spacing/2 & grid_mat(:,:,1)<=point(1)+spacing/2 ...
                            & grid_mat(:,:,2)>point(2)-spacing/2 & grid_mat(:,:,2)<=point(2)+spacing/2);
        [xTS,yTS] = ind2sub(size(tnorm),indTS);
        TSlevel(ii) = mean(mean(tnorm(xTS-1:xTS+1,yTS-1:yTS+1)));
    end
    axis off


    hold on
    subplot('Position',[0.8 0.1 0.1 0.8])
    axis tight
    caxis([tmin tmax]), axis off
%     if isempty(hc)
    hc = colorbar('West');
    set(hc,'Fontsize',18)
%     end

    if length(movieData.channels_)>1
        paxImage=movieData.channels_(2).loadImage(ii);
    %     paxImageCropped = paxImage(indULy+spacing*band:indBRy-spacing*band,indULx+spacing*band:indBRx-spacing*band);
        paxImageCropped = paxImage(grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY,grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX);
        %Scale bar
        paxImageCropped(15:16,10:10+round(2000/movieData.pixelSize_))=max(max(paxImageCropped));
        imwrite(paxImageCropped, strcat(paxPath,'/paxCroppedTif',num2str(ii,iiformat),'.tif'));
    end

    % saving
    I = getframe(h1);
    imwrite(I.cdata, strcat(tifPath,'/stressMagTif',num2str(ii,iiformat),'.tif'));

%         hgexport(h1,strcat(tifPath,'/stressMagTif',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
    hgsave(h1,strcat(figPath,'/stressMagFig',num2str(ii,iiformat)),'-v7.3')

%     set(h1,'Renderer','zbuffer')
%     saveas(h1,strcat(epsPath,'/stressMagEps',num2str(ii,iiformat),'.eps'),'eps') 
    print(h1,strcat(epsPath,'/stressMagEps',num2str(ii,iiformat),'.eps'),'-depsc2')
    hold off
%     delete(hs)
    delete(hq)
end
if pointTF
    t = (0:nFrames-1)*movieData.timeInterval_;
    figure,[AX,H1,H2] = plotyy(t,TSlevel,t,paxLevel,'plot');
    set(get(AX(1),'Ylabel'),'String','Traction Stress (Pa)') 
    set(get(AX(2),'Ylabel'),'String','Paxillin Fluorescence Intensity (A.U)') 
    xlabel('Time (sec)') 
    set(H1,'LineStyle','--')
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