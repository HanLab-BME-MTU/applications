function [pos_fltr,vec_fltr,pos_out,vec_out,id_in,id_out]=filterVectorOutliers(pos,vec,numStd,boxSizeLocFac,boxSizeGlbFac,doPlot)
% for i=1:5:75; figure(i); filterVectorOutliers(displField(i).pos,displField(i).vec,14,8,2,1); end
% tripple sort-out:
% for i=1:10:200
%     figure(1000+i);[pos_fltr,vec_fltr,pos_out,vec_out,id_in,id_out]=filterVectorOutliers(displField(i).pos,displField(i).vec,18,10,6,1);
%     figure(2000+i);[pos_fltr_2,vec_fltr_2]=filterVectorOutliers(pos_fltr,vec_fltr,18,10,6,1);
%     figure(i);[pos_fltr_3,vec_fltr_3]=filterVectorOutliers(pos_fltr_2,vec_fltr_2,18,10,6,1); 
% end

if nargin < 3 || isempty(numStd)
    % Set a threshold based on the global std:
    numStd=10;
end
if nargin < 4 || isempty(boxSizeLocFac)
    % Set the local box size (factor times the gridSize):
    boxSizeLocFac=6;
end
if nargin < 5 || isempty(boxSizeGlbFac)
    % Set the box size for the averaging (factor times the gridSize):
    boxSizeGlbFac=2;
end
if nargin < 6 || isempty(doPlot)
    doPlot=0;
end

xmin=min(pos(:,1));
xmax=max(pos(:,1));
ymin=min(pos(:,2));
ymax=max(pos(:,2));

ind = sub2ind([ymax xmax],pos(:,2),pos(:,1));

% Apply local average filter to x,y-comp. The support is chosen such that
% numPoints ~ (2+1)^2~10 (larger supports penalizes vectors 
% in the center of force nodes too much):
gridSize=ceil(sqrt((xmax-xmin)*(ymax-ymin)/length(vec)));
h = fspecial('average', round(boxSizeGlbFac*gridSize));

% interpolate the force field on all pixel positions:
ux_intp = TriScatteredInterp(pos(:,1),pos(:,2),vec(:,1),'nearest');
uy_intp = TriScatteredInterp(pos(:,1),pos(:,2),vec(:,2),'nearest');

[xgrid,ygrid] = meshgrid(1:xmax,1:ymax);
mat_x= ux_intp(xgrid,ygrid);
mat_y= uy_intp(xgrid,ygrid);

mat_x_smth = imfilter(mat_x,h);
mat_y_smth = imfilter(mat_y,h);

% The filtered values are:
vec_smth(:,1)=mat_x_smth(ind);
vec_smth(:,2)=mat_y_smth(ind);

% The deviations (noise) from the mean field are:
vec_dev=(vec-vec_smth);

% if doPlot==1
%     h=figure(3);
%     % Determine the outliers based on the global std:
%     [id_out_me,id_in_me,~] = idTwoRandVarOutlier(vec_dev(:,1),vec_dev(:,2),'threshold',numStd,'figH',h);
% else
[id_out_me,id_in_me,~] = idTwoRandVarOutlier(vec_dev(:,1),vec_dev(:,2),'threshold',numStd);
% end

% Determine the outliers based on local std:
id_out_lin = idVectorOutlier(pos,vec,boxSizeLocFac*gridSize,'lenThreshold',1,'dirThreshold',1,'pplThreshold',1);

if doPlot==1
    % my (global) results:
    pos_fltr_me=pos(id_in_me,:);
    vec_fltr_me=vec(id_in_me,:);
    pos_out_me =pos(id_out_me,:);
    vec_out_me =vec(id_out_me,:);

    % Lin (local) results:
    pos_out_lin =pos(id_out_lin,:)+1;
    vec_out_lin =vec(id_out_lin,:);

    maxForcePlot=0.3/gridSize*max(sqrt(vec(:,1).^2+vec(:,2).^2));
    %imagesc(mat_x_smth)    
    quiver(pos_fltr_me(:,1),pos_fltr_me(:,2),vec_fltr_me(:,1)/maxForcePlot,vec_fltr_me(:,2)/maxForcePlot,0,'k')
    hold on
    quiver(pos_out_me(:,1) ,pos_out_me(:,2) ,vec_out_me(:,1)/maxForcePlot ,vec_out_me(:,2)/maxForcePlot ,0,'r')
    quiver(pos_out_lin(:,1) ,pos_out_lin(:,2) ,vec_out_lin(:,1)/maxForcePlot ,vec_out_lin(:,2)/maxForcePlot ,0,'g')
    hold off

%     figure(2)
%     imagesc(mat_y_smth)
%     hold on
%     quiver(pos_fltr_me(:,1),pos_fltr_me(:,2),vec_fltr_me(:,1)/maxForcePlot,vec_fltr_me(:,2)/maxForcePlot,0,'k')
%     quiver(pos_out_me(:,1) ,pos_out_me(:,2) ,vec_out_me(:,1)/maxForcePlot ,vec_out_me(:,2)/maxForcePlot ,0,'r')
%     quiver(pos_out_lin(:,1) ,pos_out_lin(:,2) ,vec_out_lin(:,1)/maxForcePlot ,vec_out_lin(:,2)/maxForcePlot ,0,'g')
%     hold off
end

% combine the two results:
id_out=sort(union(id_out_me, id_out_lin));
id_out=id_out(:);
id_in = setdiff(id_in_me,id_out_lin);
id_in=id_in(:);

% prepare the output
pos_fltr=pos(id_in,:);
vec_fltr=vec(id_in,:);
pos_out =pos(id_out,:);
vec_out =vec(id_out,:);




% few lines of helpful code:
% for the use of a gaussian filter:
%**************************************************************************
% sigma = 15;
% h = fspecial('gaussian', 2*(5*sigma), sigma);
% 
% mat_x=0*zeros(ymax,xmax);
% mat_y=0*zeros(ymax,xmax);
% 
% mat_x(ind) = vec(:,1);
% mat_y(ind) = vec(:,2);
% 
% mat_x_smth = imfilter(mat_x,h);
% mat_y_smth = imfilter(mat_y,h);
% 
% % The filtered values are:
% vec_smth(:,1)=mat_x_smth(ind);
% vec_smth(:,2)=mat_y_smth(ind);
% 
% % This is rescaling is kind of simple minded, there has to be found a
% % better way:
% scale=median(sqrt(sum(vec.^2,2)))/median(sqrt(sum(vec_smth.^2,2)));
% mat_x_smth = scale*mat_x_smth;
% mat_y_smth = scale*mat_y_smth;
%**************************************************************************
% To estimate the local std:
% method='none';
% if strcmp(method,'local')
%     mat_x=NaN*zeros(ymax,xmax);
%     mat_y=NaN*zeros(ymax,xmax);
%     mat_x(ind) = vec(:,1);
%     mat_y(ind) = vec(:,2);
%     boxL=5*gridSize;
%     for j=1:length(vec)
%         xRange=max(1,vec(j,1)-boxL):min(xmax,vec(j,1)+boxL);
%         yRange=max(1,vec(j,2)-boxL):min(ymax,vec(j,2)+boxL);
%         
%         currMat_x=mat_x(yRange,xRange);
%         currMat_y=mat_y(yRange,xRange);
%         
%         currMat_x_smth=mat_x(yRange,xRange);
%         currMat_y_smth=mat_y(yRange,xRange);
%     end
% end
%**************************************************************************
% The function:
% idTwoRandVarOutlier(vec_dev(:,1),vec_dev(:,2),'threshold',numStd,'figH',h);
% is similar to these few lines of codes which provide the
% same results most of the time (since std_x ~ std_y because of rotational
% symmetry):
%  threshold=numStd*mean(std(vec_dev));
%  checkVec=(sqrt(sum(vec_dev.^2,2))>threshold);
%  id_out=find( checkVec);
%  id_in =find(~checkVec);