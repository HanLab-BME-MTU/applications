function tMap = generateHeatmapShifted(forceField,displField,band,imSize)

for ii=1:numel(forceField)
    reg_grid1=createRegGridFromDisplField(forceField,1,0); %2=2 times fine interpolation
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
    
% account for if displField contains more than one frame
[reg_grid,~,~,~]=createRegGridFromDisplField(displField(1),2); %2=2 times fine interpolation
for ii=1:numel(displField)
    [grid_mat,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
    grid_spacingX = grid_mat(1,2,1)-grid_mat(1,1,1);
    grid_spacingY = grid_mat(2,1,2)-grid_mat(1,1,2);
    imSizeX = grid_mat(end,end,1)-grid_mat(1,1,1)+grid_spacingX;
    imSizeY = grid_mat(end,end,2)-grid_mat(1,1,2)+grid_spacingY;
    if nargin<6
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
end

    [grid_mat,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
    pos = [reshape(grid_mat(:,:,1),[],1) reshape(grid_mat(:,:,2),[],1)]; %dense
    disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    [~,if_mat,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid);
    force_vec = [reshape(if_mat(:,:,1),[],1) reshape(if_mat(:,:,2),[],1)]; 

    [~,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],grid_mat); %1:cluster size
    tMap = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;
end

