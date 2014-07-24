function [] = cropMovieTFMPackage( pathForTheMovieDataFile )
%cropMovieTFMPackage crops movie images and displacement field for finer
%force reconstruction.
% input:    pathForTheMovieDataFile:    path to the movieData file
% output:   cropped displacement stored in pathForTheMovieDataFile/correctedDisplacementField
%           original displacement stored in pathForTheMovieDataFile/correctedDisplacementField
%           cropped images stored at a folder where their orinials were
%           original channels are moved to one folder inside its own folder
% Sangyoon Han April 2013
    if nargin==1
        band = 4;
    end
    % Load the MovieData
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
    movieData = MovieData.load(movieDataPath);
    % Get whole frame number
    nFrames = movieData.nFrames_;
    % Load the displField
    iForceFieldProc = 3;
    displFieldProc=movieData.processes_{iForceFieldProc};
    displField=displFieldProc.loadChannelOutput;
    % Set up the output file path
    outputFilePath = [pathForTheMovieDataFile filesep 'Heatmaps'];
    epsPath = [outputFilePath filesep 'eps'];
    tifPath = [outputFilePath filesep 'tifs'];
    figPath = [outputFilePath filesep 'figs'];
    if ~exist(epsPath,'dir')
        mkdir(epsPath);
        mkdir(tifPath);
        mkdir(figPath);
    end
    
    %Find the maximum force.
    tmax = 0;
    tmin = 100000;
    [reg_grid1,~,~,~]=createRegGridFromDisplField(displField,1); %2=2 times fine interpolation

    % band width for cutting border
%     band=4;

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

        tmax = max(tmax,nanmax(fnorm_vec));
        tmin = min(tmin,nanmin(fnorm_vec));
    end
    tmax = 286;
%     tmin = tmin-0.1;
%     tmax=tmax/5;
%     LeftUpperCorner(1:2) = [min(displField(1).pos(:,1)), min(displField(1).pos(:,2))];
%     RightLowerCorner(1:2) = [max(displField(1).pos(:,1)), max(displField(1).pos(:,2))];

    [reg_grid,~,~,~]=createRegGridFromDisplField(displField,3); %2=2 times fine interpolation

    h1 = figure;
    hold off
    set(h1, 'Position', [100 100 movieData.imSize_(2)*10/9 movieData.imSize_(1)])
    hc = []; %handle for colorbar
    iiformat = ['%.' '3' 'd'];

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
        grid_mat(end-band:end,:,:)=[];
        grid_mat(:,end-band:end,:)=[];
        grid_mat(1:1+band,:,:)=[];
        grid_mat(:,1:1+band,:)=[];


        % drawing
        subplot('Position',[0 0 0.9 1])
%         hs = surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm,'FaceColor','interp',...
%             'EdgeColor','none', 'FaceLighting','gouraud');%, 'FaceLighting','phong');
%         zlim([tmin tmax]), view(0,90)
        hs = pcolor(grid_mat(:,:,1), grid_mat(:,:,2), tnorm);%,[tmin tmax]);
        colormap jet;
        shading interp
        caxis([tmin tmax])
        set(gca, 'DataAspectRatio', [1,1,1],'Ydir','reverse');
%         daspect([5 5 1])

%         xlim([LeftUpperCorner(1) RightLowerCorner(1)])
%         ylim([LeftUpperCorner(2) RightLowerCorner(2)])
        axis tight
        axis off
        

        hold on
        subplot('Position',[0.9 0.1 0.1 0.8])
        axis tight
        caxis([tmin tmax]), axis off
        if isempty(hc)
            hc = colorbar('West');
        end
        
        % saving

        I = getframe(h1);
        imwrite(I.cdata, strcat(tifPath,'/stressMagTif',num2str(ii,iiformat),'.tif'));

%         hgexport(h1,strcat(tifPath,'/stressMagTif',num2str(ii,iiformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/stressMagFig',num2str(ii,iiformat)),'-v7.3')

%         print(h1,strcat(epsPath,'/stressMagEps',num2str(ii,iiformat),'.eps'),'-depsc')
        hold off
        delete(hs)
    end

return;
% to run the function:
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret-controlcell5/TFM/cell 5/c647_im',6);
