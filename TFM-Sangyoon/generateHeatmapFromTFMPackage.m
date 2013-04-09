function [] = generateHeatmapFromTFMPackage( pathForTheMovieDataFile,band )
%generateHeatmapFromTFMPackage generates heatmap from forcefield stored in
%movieData.
% input:    pathForTheMovieDataFile:    path to the movieData file
%           band:                       band width for cutting border
%           (default=4)
% output:   images of heatmap stored in pathForTheMovieDataFile/heatmap
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
    % Load the forcefield
    iForceFieldProc = 4;
    forceFieldProc=movieData.processes_{iForceFieldProc};
    forceField=forceFieldProc.loadChannelOutput;
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

        tmax = max(tmax,max(fnorm_vec));
        tmin = min(tmin,min(fnorm_vec));
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
    TSlevel = zeros(nFrames,1);

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

%         xlim([LeftUpperCorner(1) RightLowerCorner(1)])
%         ylim([LeftUpperCorner(2) RightLowerCorner(2)])
        axis tight
        % pick a point
        if ii==1
            point = ginput(1);
        end
        % point = [357.7008  319.1465]
        % grid_mat's sub for point
        spacing = grid_mat(2,1,1)-grid_mat(1,1,1);
        indTS = find(grid_mat(:,:,1)>point(1)-spacing/2 & grid_mat(:,:,1)<=point(1)+spacing/2 ...
                            & grid_mat(:,:,2)>point(2)-spacing/2 & grid_mat(:,:,2)<=point(2)+spacing/2);
        [xTS,yTS] = ind2sub(size(tnorm),indTS);
        TSlevel(ii) = mean(mean(tnorm(xTS-1:xTS+1,yTS-1:yTS+1)));
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
    figure,plot(0:nFrames-1,TSlevel,'r')

return;
% to run the function:
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret/TFM/cell 5/c647_im',6);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret/TFM/cell3/TFM',40);
generateHeatmapFromTFMPackage('/files/.retain-snapshots.d7d-w0d/LCCB/shared/X-change/forSangyoon/fromYoubean/130402 data/130110 Cell4',24);

% desktop version
generateHeatmapFromTFMPackage('/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/IntraVsExtraForce/Margaret/TFM/cell 5/c647_im',6);
