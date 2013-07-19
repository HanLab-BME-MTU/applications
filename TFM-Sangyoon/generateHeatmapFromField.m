function [umMap,XI,YI]=generateHeatmapFromField(displField,dataPath,ummax,cmapmode,w,h)
if nargin <4
    cmapmode = 'jet';
end
[reg_grid,~,~,~]=createRegGridFromDisplField(displField,2); %2=2 times fine interpolation
[grid_mat,iu_mat,~,~] = interp_vec2grid(displField.pos, displField.vec,[],reg_grid);
grid_spacingX = grid_mat(1,2,1)-grid_mat(1,1,1);
grid_spacingY = grid_mat(2,1,2)-grid_mat(1,1,2);
imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1)+grid_spacingX;
imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2)+grid_spacingY;
if nargin<5
    w = imSizeX;
    h = imSizeY;
end
centerX = ((grid_mat(end,end,1)+grid_mat(1,1,1))/2);
centerY = ((grid_mat(end,end,2)+grid_mat(1,1,2))/2);
% [XI,YI] = meshgrid(grid_mat(1,1,1):grid_mat(1,1,1)+imSizeX,grid_mat(1,1,2):grid_mat(1,1,2)+imSizeY);
xmin = centerX-w/2;
xmax = centerX+w/2;
ymin = centerY-h/2;
ymax = centerY+h/2;
[XI,YI] = meshgrid(xmin:xmax,ymin:ymax);

umnorm = (iu_mat(:,:,1).^2 + iu_mat(:,:,2).^2).^0.5;
umMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),umnorm,XI,YI,'cubic');
if nargin >=2
    ummin = min(min(umnorm));
    if nargin <3
        ummax = max(max(umnorm));
    end
    
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
        color3 = [0 1 0]; color4 = [0 0 1];
        color5 = [1 1 1]; color6 = [1 1 0];
        color7 = [0 1 1]; color8 = [1 0 1];    
        uDefinedCool = usercolormap(color1,color4,color7, color3,color6,color2);
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
    dispScale=0.25*ummax; %max(sqrt(displField.vec(:,1).^2+displField.vec(:,2).^2));
    
    Npoints = length(displField.pos(:,1));
    inIdx = false(Npoints,1);

    for ii=1:Npoints
        if displField.pos(ii,1)>=xmin && displField.pos(ii,1)<=xmax ...
                && displField.pos(ii,2)>=ymin && displField.pos(ii,2)<=ymax
            inIdx(ii) = true;
        end
    end
    
    quiver(displField.pos(inIdx,1)-xmin,displField.pos(inIdx,2)-ymin, displField.vec(inIdx,1)./dispScale,displField.vec(inIdx,2)./dispScale,0,'Color',[75/255 0/255 130/255],'LineWidth',0.5);
    
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
    hgsave(h3,strcat(figPath,'/displFieldMagFig'),'-v7.3')
    print(h3,strcat(epsPath,'/displFieldMagEps.eps'),'-depsc2')
end
end


