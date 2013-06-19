function [umMap,XI,YI]=generateHeatmapFromField(displField,dataPath,ummax)
    [reg_grid,~,~,~]=createRegGridFromDisplField(displField,2); %2=2 times fine interpolation
    [grid_mat,iu_mat,~,~] = interp_vec2grid(displField.pos, displField.vec,[],reg_grid);
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1);
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2);
    [XI,YI]=meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);

    umnorm = (iu_mat(:,:,1).^2 + iu_mat(:,:,2).^2).^0.5;
    umMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),umnorm,XI,YI,'cubic');
    if nargin >=2
        ummin = min(min(umnorm));
        if nargin <3
            ummax = max(max(umnorm));
        end

        h3=figure('color','w');
        set(h3, 'Position', [100 100 imSizeX*1.25 imSizeY])
        subplot('Position',[0 0 0.8 1])
        imshow(umMap,[ummin ummax]), colormap jet;
        %quiver plot
        hold on
        dispScale=0.25*ummax; %max(sqrt(displField.vec(:,1).^2+displField.vec(:,2).^2));

        quiver(displField.pos(:,1)-grid_mat(1,1,1),displField.pos(:,2)-grid_mat(1,1,2), displField.vec(:,1)./dispScale,displField.vec(:,2)./dispScale,0,'Color',[75/255 0/255 130/255]);

        subplot('Position',[0.8 0.1 0.1 0.8])
        axis tight
        caxis([ummin ummax]), axis off
        hc = colorbar('West');
        set(hc,'Fontsize',12)

        % saving
        % Set up the output file path
        outputFilePath = [dataPath filesep 'displacementField measured'];
        tifPath = [outputFilePath filesep 'tifs'];
        figPath = [outputFilePath filesep 'figs'];
        epsPath = [outputFilePath filesep 'eps'];
        if ~exist(tifPath,'dir') || ~exist(epsPath,'dir')
            mkdir(tifPath);
            mkdir(figPath);
            mkdir(epsPath);
        end

        I = getframe(h3);
        imwrite(I.cdata, strcat(tifPath,'/displFieldMagTif','.tif'));
        hgsave(h3,strcat(figPath,'/displFieldMagFig','-v7.3'))
        print(h3,strcat(epsPath,'/displFieldMagEps.eps'),'-depsc2')
    end
end
    

