function [umMap,XI,YI]=generateHeatmapFromFieldBare(Field,dataPath,band,ummax,cmapmode,w,h,plotQuiver)
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
    umMapX = griddata(grid_mat(:,:,1),grid_mat(:,:,2),iu_mat(:,:,1),XI,YI,'cubic');
    umMapY = griddata(grid_mat(:,:,1),grid_mat(:,:,2),iu_mat(:,:,2),XI,YI,'cubic');

    % saving
    % Set up the output file path
    if ~isempty(dataPath)
        outputFilePath = [dataPath filesep 'heatMap'];
        if ~exist(outputFilePath,'dir')
            mkdir(outputFilePath);
        end

        imwrite(uint16(umMap), strcat(outputFilePath,'/fieldMag',num2str(k),'.tif'));
        imwrite(uint16(umMapX), strcat(outputFilePath,'/fieldMagX',num2str(k),'.tif'));
        imwrite(uint16(umMapY), strcat(outputFilePath,'/fieldMagY',num2str(k),'.tif'));
    end
end
end


