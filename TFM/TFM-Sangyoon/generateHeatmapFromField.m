function [umMap,XI,YI]=generateHeatmapFromField(Field,dataPath,band,ummax,cmapmode,w,h,plotQuiver)
if nargin <5 || isempty(cmapmode)
    cmapmode = 'jet';
end
if nargin <8
    plotQuiver = true;
end
ummin = 1e20;
if nargin <4 || isempty(ummax)
    temp_ummax = 0;
end
if nargin <3
    band = 0;
end
numNonEmpty = sum(arrayfun(@(x) ~isempty(x.vec),Field));
for k=1:numNonEmpty
    maxMag = (Field(k).vec(:,1).^2+Field(k).vec(:,2).^2).^0.5;
    ummin = min(ummin,min(maxMag));
    if nargin < 4 || isempty(ummax)
        temp_ummax = max(temp_ummax, max(maxMag));
    end
end
if nargin < 4 || isempty(ummax)
    ummax = temp_ummax;
end

% account for if displField contains more than one frame
[reg_grid,~,~,~]=createRegGridFromDisplField(Field(1),2); %2=2 times fine interpolation
for k=1:numNonEmpty
    [grid_mat,iu_mat,~,~] = interp_vec2grid(Field(k).pos, Field(k).vec,[],reg_grid);
    grid_spacingX = grid_mat(1,2,1)-grid_mat(1,1,1);
    grid_spacingY = grid_mat(2,1,2)-grid_mat(1,1,2);
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1)+grid_spacingX;
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2)+grid_spacingY;
    if nargin<6 || isempty(w) || isempty(h)
        w = imSizeX;
        h = imSizeY;
    end
    centerX = ((grid_mat(end,end,1)+grid_mat(1,1,1))/2);
    centerY = ((grid_mat(end,end,2)+grid_mat(1,1,2))/2);
    % [XI,YI] = meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
    xmin = centerX-w/2+band;
    xmax = centerX+w/2-band;
    ymin = centerY-h/2+band;
    ymax = centerY+h/2-band;
    [XI,YI] = meshgrid(xmin:xmax,ymin:ymax);

    umnorm = (iu_mat(:,:,1).^2 + iu_mat(:,:,2).^2).^0.5;
    umMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),umnorm,XI,YI,'cubic');
    if nargin >=2
        h3=figure('color','w');
        set(h3, 'Position', [100 100 w*1.25 h])
        subplot('Position',[0 0 0.8 1])
        imshow(umMap,[ummin ummax])
        if strcmp(cmapmode,'uDefinedCool')
            color1 = [0 0 0]; color2 = [1 0 0];
            color3 = [0 1 0]; color4 = [0 0 1];
            color5 = [1 1 1]; color6 = [1 1 0];
            color7 = [0 1 1]; color8 = [1 0 1];    
            uDefinedCool = usercolormap(color1,color4,color7,color5);
            colormap(uDefinedCool);
        elseif strcmp(cmapmode,'uDefinedJet')
            color1 = [0 0 0]; color2 = [1 0 0];
            color3 = [122/255 179/255 23/255]; color4 = [0 0 1];
            color5 = [1 1 1]; color6 = [1 1 0];
            color7 = [0 1 1]; color8 = [252/255 145/255 58/255];    
            uDefinedCool = usercolormap(color1,color4,color7, color3,color6,color8,color2);
            colormap(uDefinedCool);
        elseif strcmp(cmapmode,'uDefinedRYG')
            color1 = [0 0 0]; color2 = [1 0 0];
            color3 = [0 1 0]; color4 = [0 0 1];
            color5 = [1 1 1]; color6 = [1 1 0];
            color7 = [0 1 1]; color8 = [1 0 1];   
            color9 = [49/255 0 98/255];
            uDefinedCool = usercolormap(color9,color2,color6, color3);
            colormap(uDefinedCool);
        elseif strcmp(cmapmode,'uDefinedYGR')
            color1 = [0 0 0]; color2 = [1 0 0];
            color3 = [0 1 0]; color4 = [0 0 1];
            color5 = [1 1 1]; color6 = [1 1 0];
            color7 = [0 1 1]; color8 = [1 0 1];    
            uDefinedCool = usercolormap(color5,color6,color3,color2);
            colormap(uDefinedCool);
        else
            colormap(cmapmode);
        end
        %quiver plot
        hold on
        dispScale=0.05*ummax; %max(sqrt(displField.vec(:,1).^2+displField.vec(:,2).^2));

        Npoints = length(Field(k).pos(:,1));
        inIdx = false(Npoints,1);

        for ii=1:Npoints
            if Field(k).pos(ii,1)>=xmin && Field(k).pos(ii,1)<=xmax ...
                    && Field(k).pos(ii,2)>=ymin && Field(k).pos(ii,2)<=ymax
                inIdx(ii) = true;
            end
        end

        if plotQuiver
            quiver(Field(k).pos(inIdx,1)-xmin,Field(k).pos(inIdx,2)-ymin, Field(k).vec(inIdx,1)./dispScale,Field(k).vec(inIdx,2)./dispScale,0,'Color',[75/255 0/255 130/255],'LineWidth',0.5);
        end

        subplot('Position',[0.8 0.1 0.1 0.8])
        axis tight
        caxis([ummin ummax]), axis off
        hc = colorbar('West');
        set(hc,'Fontsize',12)

        % saving
        % Set up the output file path
        if ~isempty(dataPath)
            outputFilePath = [dataPath filesep 'heatMap'];
            tifPath = [outputFilePath filesep 'tifs'];
            figPath = [outputFilePath filesep 'figs'];
            epsPath = [outputFilePath filesep 'eps'];
            if ~exist(tifPath,'dir') || ~exist(epsPath,'dir')
                mkdir(tifPath);
                mkdir(figPath);
                mkdir(epsPath);
            end

            I = getframe(h3);
            imwrite(I.cdata, strcat(tifPath,'/displFieldMagTif',num2str(k),'.tif'));
            hgsave(h3,strcat(figPath,'/displFieldMagFig',num2str(k)),'-v7.3')
            print(h3,strcat(epsPath,'/displFieldMagEps',num2str(k),'.eps'),'-depsc2')
            close(h3)
        end
    end
end
end


