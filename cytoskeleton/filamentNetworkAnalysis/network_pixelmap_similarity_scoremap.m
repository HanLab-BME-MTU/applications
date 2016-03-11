function [similarity_scoremap, difference_map] = ...
    network_pixelmap_similarity_scoremap(VIF_orient_map,MT_orient_map, radius)
% Function for the network comparison with pixel based input( no filament info) version

%% for debugging with traction force orientation
% VIF_orient_map = -D;
% VIF_orient_map = rand(size(current_seg_orientation)/25)*pi-pi/2;
% VIF_orient_map = imresize(VIF_orient_map,size(current_seg_orientation));
% VIF_orient_map(isnan(current_seg_orientation)>0)=nan;
% MT_orient_map = current_seg_orientation;
% radius=30;


%%
img_size = size(VIF_orient_map);

distance_map_1_2 = nan(img_size);
distance_map_2_1 = nan(img_size);
angle_map_1_2 = nan(img_size);
angle_map_2_1 = nan(img_size);
angle_map_1_2_A = nan(img_size);
angle_map_2_1_A = nan(img_size);
angle_map_1_2_B = nan(img_size);
angle_map_2_1_B = nan(img_size);

angle_map_II_1_2_A= nan(img_size);
angle_map_II_1_2_B= nan(img_size);
angle_map_II_2_1_A= nan(img_size);
angle_map_II_2_1_B= nan(img_size);

orient_map_1_2_A = nan(img_size);
orient_map_1_2_B = nan(img_size);
orient_map_2_1_A = nan(img_size);
orient_map_2_1_B = nan(img_size);
orient_map_1_2 = nan(img_size);
orient_map_2_1 = nan(img_size);
%%

VIF_current_seg = (isnan(VIF_orient_map)==0);
MT_current_seg = (isnan(MT_orient_map)==0);


[VIF_YY,VIF_XX] = find(isnan(VIF_orient_map)==0);
VIF_OO = VIF_orient_map(isnan(VIF_orient_map)==0);
VIF_II = ones(size(VIF_YY));


[MT_YY, MT_XX] = find(isnan(MT_orient_map)==0);
MT_OO = MT_orient_map(isnan(MT_orient_map)==0);
MT_II = ones(size(MT_YY));




% VIF as channel 1, MT as channel 2, use 1 and 2 for shortness
X1=VIF_XX;
Y1=VIF_YY;

X2=MT_XX;
Y2=MT_YY;

O1 = VIF_OO;
O2 = MT_OO;


% find in the pool of the points in the 2nd channel, the closest point of
% the points in the 1st channel, so query is 1, pool is 2.
% [idx_cell, dist_cell] = KDTreeBallQuery([Y2 X2],[Y1 X1],radius);

% matlab kdtree built in code
[idx_cell, dist_cell] = rangesearch([Y2 X2],[Y1 X1],radius,'nsmethod','kdtree');


for iQ = 1 : length(Y1)
    idx_this = idx_cell{iQ}';
    dist_this = dist_cell{iQ}';
    
    [sort_dist,sort_IX] = sort(dist_this);
    
    if(length(sort_dist)>0)
        sort_idx = idx_this(sort_IX);
        
        distance_map_1_2(sub2ind(img_size,Y1(iQ),X1(iQ)))=sort_dist(1);
        
        %  angle is the direction difference
        angle_map_1_2_A(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
            VIF_OO(iQ) - MT_OO(sort_idx(1));
        
        % orient is the direction 
        orient_map_1_2_A(sub2ind(img_size,Y1(iQ),X1(iQ))) = MT_OO(sort_idx(1));
        
        angle_map_II_1_2_A(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
            MT_II(sort_idx(1));
        
        % if found two point, same point, could be crossing, get the other one as well.
        if(length(sort_dist)>=2 && length(sort_idx)>=2)
            if sort_dist(1)==sort_dist(2) && Y2(sort_idx(1))==Y2(sort_idx(2))...
                    && X2(sort_idx(1))==X2(sort_idx(2))
                angle_map_1_2_B(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
                    VIF_OO(iQ) - MT_OO(sort_idx(2));
                angle_map_II_1_2_B(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
                    MT_II(sort_idx(2));
                orient_map_1_2_B(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
                    MT_OO(sort_idx(2));
        
            end
        end
    end
end

angle_map_1_2_A(angle_map_1_2_A>pi/2) = angle_map_1_2_A(angle_map_1_2_A>pi/2) -pi;
angle_map_1_2_A(angle_map_1_2_A<-pi/2) = angle_map_1_2_A(angle_map_1_2_A<-pi/2) +pi;

angle_map_1_2_A(angle_map_1_2_A>pi/2) = angle_map_1_2_A(angle_map_1_2_A>pi/2) -pi;
angle_map_1_2_A(angle_map_1_2_A<-pi/2) = angle_map_1_2_A(angle_map_1_2_A<-pi/2) +pi;

angle_map_1_2_B(angle_map_1_2_B>pi/2) = angle_map_1_2_B(angle_map_1_2_B>pi/2) -pi;
angle_map_1_2_B(angle_map_1_2_B<-pi/2) = angle_map_1_2_B(angle_map_1_2_B<-pi/2) +pi;

angle_map_1_2_B(angle_map_1_2_B>pi/2) = angle_map_1_2_B(angle_map_1_2_B>pi/2) -pi;
angle_map_1_2_B(angle_map_1_2_B<-pi/2) = angle_map_1_2_B(angle_map_1_2_B<-pi/2) +pi;


angle_map_1_2 = nan(size(angle_map_1_2_A,1),size(angle_map_1_2_A,2));
map_II_1_2 = nan(size(angle_map_1_2_A,1),size(angle_map_1_2_A,2));
 

for i = 1 : size(angle_map_1_2_A, 1)
    for j = 1 : size(angle_map_1_2_A, 2)
        % if there is data
        if(~isnan(angle_map_1_2_A(i,j)))
            % if there is only one closest point
            if(isnan(angle_map_1_2_B(i,j)))
                map_II_1_2(i,j)=angle_map_II_1_2_A(i,j);
                angle_map_1_2(i,j)= abs(angle_map_1_2_A(i,j));
                
                orient_map_1_2(i,j) = orient_map_1_2_A(i,j);
            else
                % if it is closest to a crossing(2 points), pick the one
                % with smaller abs value
                if(abs(angle_map_1_2_A(i,j))>abs(angle_map_1_2_B(i,j)))
                    map_II_1_2(i,j)=angle_map_II_1_2_B(i,j);
                    angle_map_1_2(i,j)= abs(angle_map_1_2_B(i,j));  
                    orient_map_1_2(i,j) = orient_map_1_2_B(i,j);

                else
                    map_II_1_2(i,j)=angle_map_II_1_2_A(i,j);
                    angle_map_1_2(i,j)= abs(angle_map_1_2_A(i,j));   
                    orient_map_1_2(i,j) = orient_map_1_2_A(i,j);

                end % end picking
            end % end if only one point
        end % end if there is data
    end % end j
end % end i
    

% % if the dist is farther than radius, gate it, so this index is for 1
% ind_gate = find(dist<=radius);
% 
% % then update the index for founded closest points
% idx_close = idx(ind_gate);
% 
% % screen out the gated ones
% X1_close = X1(ind_gate);
% Y1_close = Y1(ind_gate);
% X2_close = X2(idx_close);
% Y2_close = Y2(idx_close);
% 
% dist_gate = dist(ind_gate);
% 
% distance_map_1_2(sub2ind(img_size,Y1_close,X1_close))=dist_gate;
% angle_map_1_2(sub2ind(img_size,Y1_close,X1_close)) = ...
%     VIF_orientation(sub2ind(img_size,Y1_close,X1_close)) ...
%     - MT_orientation(sub2ind(img_size,Y2_close,X2_close));
% 

 %%
%then the other way around 2->1

% [idx_cell, dist_cell] = KDTreeBallQuery([Y1 X1],[Y2 X2],radius);
% matlab kdtree built in code
[idx_cell, dist_cell] = rangesearch([Y1 X1],[Y2 X2],radius,'nsmethod','kdtree');


for iQ = 1 : length(Y2)
    idx_this = idx_cell{iQ}';
    dist_this = dist_cell{iQ}';
    
    [sort_dist,sort_IX] = sort(dist_this);
    
    if(length(sort_dist)>0)
        
        sort_idx = idx_this(sort_IX);
        if(length(sort_dist)>0)
            distance_map_2_1(sub2ind(img_size,Y2(iQ),X2(iQ)))=sort_dist(1);
            angle_map_2_1_A(sub2ind(img_size,Y2(iQ),X2(iQ))) = ...
                MT_OO(iQ) - VIF_OO(sort_idx(1));
            angle_map_II_2_1_A(sub2ind(img_size,Y2(iQ),X2(iQ))) = ...
                VIF_II(sort_idx(1));
            orient_map_2_1_A(i,j) = VIF_OO(sort_idx(1));

        end
        
        if(length(sort_dist)>=2 && length(sort_idx)>=2)
            if sort_dist(1)==sort_dist(2) && Y1(sort_idx(1))==Y1(sort_idx(2))...
                    && X1(sort_idx(1))==X1(sort_idx(2))
                angle_map_2_1_B(sub2ind(img_size,Y2(iQ),X2(iQ))) = ...
                    MT_OO(iQ) - VIF_OO(sort_idx(2));
                angle_map_II_2_1_B(sub2ind(img_size,Y2(iQ),X2(iQ))) = ...
                    VIF_II(sort_idx(2));
                orient_map_2_1_B(i,j) = VIF_OO(sort_idx(2));

            end
        end
    end
end

angle_map_2_1_A(angle_map_2_1_A>pi/2) = angle_map_2_1_A(angle_map_2_1_A>pi/2) -pi;
angle_map_2_1_A(angle_map_2_1_A<-pi/2) = angle_map_2_1_A(angle_map_2_1_A<-pi/2) +pi;

angle_map_2_1_A(angle_map_2_1_A>pi/2) = angle_map_2_1_A(angle_map_2_1_A>pi/2) -pi;
angle_map_2_1_A(angle_map_2_1_A<-pi/2) = angle_map_2_1_A(angle_map_2_1_A<-pi/2) +pi;

angle_map_2_1_B(angle_map_2_1_B>pi/2) = angle_map_2_1_B(angle_map_2_1_B>pi/2) -pi;
angle_map_2_1_B(angle_map_2_1_B<-pi/2) = angle_map_2_1_B(angle_map_2_1_B<-pi/2) +pi;

angle_map_2_1_B(angle_map_2_1_B>pi/2) = angle_map_2_1_B(angle_map_2_1_B>pi/2) -pi;
angle_map_2_1_B(angle_map_2_1_B<-pi/2) = angle_map_2_1_B(angle_map_2_1_B<-pi/2) +pi;

% angle_map_2_1 = min(abs(angle_map_2_1_A),abs(angle_map_2_1_B));


map_II_2_1 = nan(size(angle_map_2_1_A,1),size(angle_map_2_1_A,2));
angle_map_2_1 = nan(size(angle_map_2_1_A,1),size(angle_map_2_1_A,2));
 
for i = 1 : size(angle_map_2_1_A, 1)
    for j = 1 : size(angle_map_2_1_A, 2)
        % if there is data
        if(~isnan(angle_map_2_1_A(i,j)))
            % if there is only one closest point
            if(isnan(angle_map_2_1_B(i,j)))
                map_II_2_1(i,j)=angle_map_II_2_1_A(i,j);
                angle_map_2_1(i,j)= abs(angle_map_2_1_A(i,j));
                  orient_map_2_1(i,j) = orient_map_2_1_A(i,j);

            else
                % if it is closest to a crossing(2 points), pick the one
                % with smaller abs value
                if(abs(angle_map_2_1_A(i,j))>abs(angle_map_2_1_B(i,j)))
                    map_II_2_1(i,j)=angle_map_II_2_1_B(i,j);
                    angle_map_2_1(i,j)= abs(angle_map_2_1_B(i,j));  
                    orient_map_2_1(i,j) = orient_map_2_1_B(i,j);
                else
                    map_II_2_1(i,j)=angle_map_II_2_1_A(i,j);
                    angle_map_2_1(i,j)= abs(angle_map_2_1_A(i,j)); 
                    orient_map_2_1(i,j) = orient_map_2_1_A(i,j);

                end % end picking
            end % end if only one point
        end % end if there is data
    end % end j
end % end i


% 
% 
% % find in the pool of the points in the 1nd channel, the closest point of
% % the points in the 2st channel, so query is 2, pool is 1 .
% [idx, dist] = KDTreeClosestPoint([Y1 X1],[Y2 X2]);
% 
% % if the dist is farther than radius, gate it, so this index is for 2
% ind_gate = find(dist<=radius);
% 
% % then update the index for founded closest points
% idx_close = idx(ind_gate);
% 
% % screen out the gated ones
% X2_close = X2(ind_gate);
% Y2_close = Y2(ind_gate);
% X1_close = X1(idx_close);
% Y1_close = Y1(idx_close);
% 
% dist_gate = dist(ind_gate);
% 
% distance_map_2_1(sub2ind(img_size,Y2_close,X2_close)) = dist_gate;
% angle_map_2_1(sub2ind(img_size,Y2_close,X2_close)) = ...
%     MT_orientation(sub2ind(img_size,Y2_close,X2_close)) ...
%     - VIF_orientation(sub2ind(img_size,Y1_close,X1_close));
% 
% 
% angle_map_1_2(angle_map_1_2>pi/2) = angle_map_1_2(angle_map_1_2>pi/2) -pi;
% angle_map_1_2(angle_map_1_2<-pi/2) = angle_map_1_2(angle_map_1_2<-pi/2) +pi;
% 
% 
% angle_map_2_1(angle_map_2_1>pi/2) = angle_map_2_1(angle_map_2_1>pi/2) -pi;
% angle_map_2_1(angle_map_2_1<-pi/2) = angle_map_2_1(angle_map_2_1<-pi/2) +pi;
% 
% 
% 
% imwrite(distance_map_1_2,[outdir,filesep,'Dis12_frame_',num2str(iFrame),'.tif']);  end;
% imwrite(distance_map_2_1,[outdir,filesep,'Dis21_frame_',num2str(iFrame),'.tif']);  end;
% imwrite(angle_map_2_1,[outdir,filesep,'Ang21_frame_',num2str(iFrame),'.tif']);  end;
% imwrite(angle_map_1_2,[outdir,filesep,'Ang12_frame_',num2str(iFrame),'.tif']);  end;
% 
whole_ROI = imdilate(VIF_current_seg,ones(radius,radius)) + imdilate(MT_current_seg,ones(radius,radius))>0;
% 
% % disgard the boundary ones
% whole_ROI(1:radius+1,:)=0;
% whole_ROI(end-radius:end,:)=0;
% whole_ROI(:,1:radius+1)=0;
% whole_ROI(:,end-radius:end)=0;

% find the points of interest
[Y,X] = find(whole_ROI>0);

score_maps_distance_1_2 = nan(img_size);
score_maps_distance_2_1 = nan(img_size);
score_maps_angle_1_2 = nan(img_size);
score_maps_angle_2_1 = nan(img_size);

[cy,cx] = find(fspecial('disk',radius)>0);


distance_map_1_2_pad = nan(img_size+2*radius);
distance_map_2_1_pad = nan(img_size+2*radius);
angle_map_1_2_pad = nan(img_size+2*radius);
angle_map_2_1_pad = nan(img_size+2*radius);


distance_map_1_2_pad(radius+1:end-radius,radius+1:end-radius) = distance_map_1_2;
distance_map_2_1_pad(radius+1:end-radius,radius+1:end-radius) = distance_map_2_1;
angle_map_1_2_pad(radius+1:end-radius,radius+1:end-radius) = angle_map_1_2;
angle_map_2_1_pad(radius+1:end-radius,radius+1:end-radius) = angle_map_2_1;

% 
% distance_map_1_2_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% distance_map_2_1_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% angle_map_1_2_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% angle_map_2_1_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% 
% 
% for s = 1 : length(cx)
%     distance_map_1_2_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = distance_map_1_2;
%     distance_map_2_1_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = distance_map_2_1;
%     angle_map_1_2_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = angle_map_1_2;
%     angle_map_2_1_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = angle_map_2_1;
% end
% 
% 
Weight_mask = fspecial('gaussian',2*radius+1,(radius*1.3)/(2)/2);
Weight_mask = Weight_mask(sub2ind([2*radius+1,2*radius+1,],cy,cx));

% for all these points, get local support
for j = 1 : length(Y)
%     j
    x = X(j);
    y = Y(j);
        
    dis_1 = distance_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));
    dis_2 = distance_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));

    ang_1 = abs(angle_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1)));
    ang_2 = abs(angle_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1)));
    
    if(~isempty(dis_1))
        score_maps_distance_1_2(y,x) = sum(dis_1(~isnan(dis_1)).*Weight_mask(~isnan(dis_1)))/sum(Weight_mask(~isnan(dis_1)));
        score_maps_angle_1_2(y,x) = sum(ang_1(~isnan(ang_1)).*Weight_mask(~isnan(ang_1)))/sum(Weight_mask(~isnan(ang_1)));
    end
    if(~isempty(dis_2))
        score_maps_distance_2_1(y,x) = sum(dis_2(~isnan(dis_2)).*Weight_mask(~isnan(dis_2)))/sum(Weight_mask(~isnan(dis_2)));       
        score_maps_angle_2_1(y,x) = sum(ang_2(~isnan(ang_2)).*Weight_mask(~isnan(ang_2)))/sum(Weight_mask(~isnan(ang_2)));
    end    
end

% fill the empty locations
score_maps_distance_2_1(isnan(score_maps_distance_2_1))=radius;
score_maps_distance_1_2(isnan(score_maps_distance_1_2))=radius;
score_maps_angle_2_1(isnan(score_maps_angle_2_1))=pi/2;
score_maps_angle_1_2(isnan(score_maps_angle_1_2))=pi/2;



% calculation of similarity score.
similarity_scoremap = exp(-(score_maps_distance_2_1+score_maps_distance_1_2).^2/(((radius*1.5)/2*sqrt(2))^2))...
    .*exp(-(abs(score_maps_angle_2_1/2)+abs(score_maps_angle_1_2/2)).^2/(1.5*(pi/3)^2));

% calculation of similarity score only consider 1->2.
similarity_scoremap_1to2 = exp(-(0+score_maps_distance_1_2*2).^2/(((radius*1.5)/2*sqrt(2))^2))...
    .*exp(-(abs(0/2)+abs(score_maps_angle_1_2/2*2)).^2/(1.5*(pi/3)^2));

% calculation of similarity score only consider 2->1.
similarity_scoremap_2to1 = exp(-(score_maps_distance_2_1*2+0).^2/(((radius*1.5)/2*sqrt(2))^2))...
    .*exp(-(abs(score_maps_angle_2_1/2*2)+abs(0/2)).^2/(1.5*(pi/3)^2));


difference_map.similarity_scoremap_proximity = exp(-(score_maps_distance_2_1+score_maps_distance_1_2).^2/(((radius*1.5)/2*sqrt(2))^2));
difference_map.similarity_scoremap_alignment = exp(-(abs(score_maps_angle_2_1/2)+abs(score_maps_angle_1_2/2)).^2/(1.5*(pi/3)^2));
difference_map.similarity_scoremap_combined = similarity_scoremap;

difference_map.similarity_scoremap_1to2 = similarity_scoremap_1to2;
difference_map.similarity_scoremap_2to1 = similarity_scoremap_2to1;
difference_map.distance_map_1_2 = distance_map_1_2;
difference_map.distance_map_2_1 = distance_map_2_1;
difference_map.angle_map_1_2 = angle_map_1_2;
difference_map.angle_map_2_1 = angle_map_2_1;
difference_map.orient_map_1_2 = orient_map_1_2;
difference_map.orient_map_2_1 = orient_map_2_1;
difference_map.score_maps_distance_1_2 = score_maps_distance_1_2;
difference_map.score_maps_distance_2_1 = score_maps_distance_2_1;
difference_map.score_maps_angle_1_2 = score_maps_angle_1_2;
difference_map.score_maps_angle_2_1 = score_maps_angle_2_1;

