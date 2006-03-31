clear
mpm_variable = load('mpm.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine if ther is MPM and M variable
% stored or only MPM
a=char(fieldnames(mpm_variable));
if size(a,1) > 1
    MPM = mpm_variable.(a(2,:));
else
    MPM = mpm_variable.(a);
end
clear mpm_variable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trace MPM
t_n =1;
track_length(1)=0;


strg = sprintf('%%.%dd',5);
backSpc = ['\b\b\b\b\b'];
    
size(MPM,1)
track_data=zeros(100000,4);
track_length=zeros(100000,1);

for i=1:size(MPM,1)
    
    fprintf(1,[strg],i); 
    
    j=1;
    while j < size(MPM,2)-1 
        if MPM(i,j) ~= 0
            j_start =j;

            if MPM(i,j+2) ~= 0
                % save start point of track
                track_data(t_n,1) = MPM(i,j);
                track_data(t_n,2) = MPM(i,j+1);

                while j < size(MPM,2) && MPM(i,j) ~= 0
                    track_length(t_n) = track_length(t_n) + 1;

                    % save endpoint of track
                    end_track_1 = MPM(i,j);
                    end_track_2 = MPM(i,j+1);
                    if  end_track_1 == 0 | end_track_2 == 0
                        disp('there is some problem');
                    end
                    j=j+2;
                end

                % save life time of track
                track_data(t_n,3) = track_length(t_n);
                % save average track velocity
                track_data(t_n,4)= sqrt((end_track_1- track_data(t_n,1))^2 +...
                    (end_track_2- track_data(t_n,2))^2)/track_length(t_n);

                t_n = t_n + 1;
                track_length(t_n) = 0;
            end
            j=j+2;   
        else
            j=j+2;
        end
    end
    j=1;
    
    fprintf(1,backSpc);
end

% remove
track_data(find(~track_data(:,1)),:) = [];

mkdir('track_analysis');

% plot life time - velocity scatter plot
h=figure;
plot(track_data(:,3),track_data(:,4),'.');
xlabel('Life time (frames)');
ylabel('Velocity (pixel/frame)');
print(h, 'track_analysis/lifeTime_velocity.tif','-dtiff');
hgsave(h, 'track_analysis/lifeTime_velocity.fig');



life_time_bins = 0.5:20.5;
velocity_bins  = 0:0.25:7;
for life = 1:20
    for vel = 1:28
       num_sp(life,vel) = sum(inpolygon(track_data(:,3),track_data(:,4),[life_time_bins(life);  life_time_bins(life+1);...
                                                                    life_time_bins(life+1);life_time_bins(life)],...
                                                                    [velocity_bins(vel);    velocity_bins(vel);...
                                                                    velocity_bins(vel+1);  velocity_bins(vel+1)]));
    end
end
life_time_bins(end) = [];
velocity_bins(end) = [];
figure,h=surface(velocity_bins, life_time_bins,num_sp);
set(h,'Marker','o');
set(h,'MarkerSize', 0.1);
set(h,'FaceColor','interp');
%set(h,'EdgeColor','none');
xlabel('Velocity (pixel/frame)');
ylabel('Life time (frames)');
view(2)
hgsave(h, 'track_analysis/lifeTime_velocity_colorplot.fig');
%figure,surface(num_sp);





% plot track life time map
% track_map   = zeros(max(track_data(:,1)),max(track_data(:,2)));
% track_map_s = zeros(max(track_data(:,1)),max(track_data(:,2)));
% track_map_l = zeros(max(track_data(:,1)),max(track_data(:,2)));
% max_length = max(track_data(:,3));
% cmap = colormap(jet(max_length));
% cmap = cmap;
% %cmap = colormap(jet(30));
% for i=1:size(track_data,1)
%      track_map(track_data(i,1),track_data(i,2),1) = cmap(track_data(i,3),1);
%      track_map(track_data(i,1),track_data(i,2),2) = cmap(track_data(i,3),2);
%      track_map(track_data(i,1),track_data(i,2),3) = cmap(track_data(i,3),3);
%      %track_map(track_data(i,1),track_data(i,2),1) = 1;
%      %track_map(track_data(i,1),track_data(i,2),2) = 1;
%      %track_map(track_data(i,1),track_data(i,2),3) = 1;
% end
% figure,imshow(track_map,[]);

threshold = 20;
track_data_short = track_data(track_data(:,3) <= threshold,:);
track_data_short = track_data_short(track_data_short(:,3) >  2,:);
track_data_long  = track_data(track_data(:,3) > threshold, :);
for i=1:size(track_data_short,1)
     track_map_short(int16(track_data_short(i,1)),int16(track_data_short(i,2)),1) = 1;
     track_map_short(int16(track_data_short(i,1)),int16(track_data_short(i,2)),2) = 0;
     track_map_short(int16(track_data_short(i,1)),int16(track_data_short(i,2)),3) = 0;
end
for i=1:size(track_data_long,1)
     track_map_long(track_data_long(i,1),track_data_long(i,2),1) = 1;
     track_map_long(track_data_long(i,1),track_data_long(i,2),2) = 0;
     track_map_long(track_data_long(i,1),track_data_long(i,2),3) = 1;
end
figure,imshow(track_map_short,[]);
imwrite(track_map_short,'track_analysis/track_map_short.tif','tif');
figure,imshow(track_map_long,[]);
imwrite(track_map_long,'track_analysis/track_map_long.tif','tif');



vel_thresh  = 2;
track_map_slow = zeros(max(track_data(:,1)),max(track_data(:,2)));
track_map_fast = zeros(max(track_data(:,1)),max(track_data(:,2)));

track_data_slow  = track_data(track_data(:,4) <= vel_thresh,:);
track_data_fast  = track_data(track_data(:,4) >  vel_thresh, :);
for i=1:size(track_data_slow,1)
     track_map_slow(int16(track_data_slow(i,1)),int16(track_data_slow(i,2)),1) = 1;
     track_map_slow(int16(track_data_slow(i,1)),int16(track_data_slow(i,2)),2) = 0;
     track_map_slow(int16(track_data_slow(i,1)),int16(track_data_slow(i,2)),3) = 0;
end
for i=1:size(track_data_fast,1)
     track_map_fast(track_data_fast(i,1),track_data_fast(i,2),1) = 1;
     track_map_fast(track_data_fast(i,1),track_data_fast(i,2),2) = 0;
     track_map_fast(track_data_fast(i,1),track_data_fast(i,2),3) = 1;
end
figure,imshow(track_map_slow,[]);
imwrite(track_map_slow,'track_analysis/track_map_slow.tif','tif');
figure,imshow(track_map_fast,[]);
imwrite(track_map_fast,'track_analysis/track_map_fast.tif','tif');


% save all data
save('track_analysis/track_data','track_data');
slow_sores = ones(size(track_data_slow));
fast_sores = ones(size(track_data_fast));
slow_scores(:,2:3) = track_data_slow(:,1:2);
slow_scores(:,4)   = track_data_slow(:,4);
save('track_analysis/slow_scores','slow_scores');
fast_scores(:,2:3) = track_data_fast(:,1:2);
fast_scores(:,4)   = track_data_fast(:,4);
save('track_analysis/fast_scores','fast_scores');


t_n
[hi x_bin] = hist(track_length,40);
figure
plot(x_bin,hi);



