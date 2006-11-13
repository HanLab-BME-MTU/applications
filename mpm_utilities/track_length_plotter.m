function track_length_plotter(MPM, result_Dir)
%mpm_variable = load('mpm.mat');

if isempty(MPM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % determine if ther is MPM and M variable
    % stored or only MPM
    a=char(fieldnames(mpm_variable));
    if size(a,1) > 1
        MPM = mpm_variable.(a(2,:));
    else
        MPM = mpm_variable.(a);
    end
    clear mpm_variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

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
                    (end_track_2- track_data(t_n,2))^2)/(track_length(t_n)-1);

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

% life_time_bins = 0.5:20.5;
% velocity_bins  = 0:0.25:7;
% for life = 1:20
%     for vel = 1:28
%        num_sp(life,vel) = sum(inpolygon(track_data(:,3),track_data(:,4),[life_time_bins(life);  life_time_bins(life+1);...
%                                                                     life_time_bins(life+1);life_time_bins(life)],...
%                                                                     [velocity_bins(vel);    velocity_bins(vel);...
%                                                                     velocity_bins(vel+1);  velocity_bins(vel+1)]));
%     end
% end
% life_time_bins(end) = [];
% velocity_bins(end) = [];
% figure,h=surface(velocity_bins, life_time_bins,num_sp);
% set(h,'Marker','o');
% set(h,'MarkerSize', 0.1);
% set(h,'FaceColor','interp');
% %set(h,'EdgeColor','none');
% xlabel('Velocity (pixel/frame)');
% ylabel('Life time (frames)');
% view(2)
% hgsave(h, [result_Dir filesep 'lifeTime_velocity_colorplot.fig']);
%figure,surface(num_sp);

track_data(find(~track_data(:,4)),:) = [];
% plot life time - velocity scatter plot

scrsz = get(0,'ScreenSize');
h=figure('Position',[10 scrsz(4)/2-70 scrsz(3)/2 scrsz(4)/2]);
plot(track_data(:,3),track_data(:,4),'.');
xlabel('Life time (frames)');
ylabel('Velocity (pixel/frame)');

print(h, [result_Dir filesep 'lifeTime_velocity.tif'],'-dtiff');
print(h, [result_Dir filesep 'lifeTime_velocity.eps'],'-depsc2','-tiff'); 
hgsave(h, [result_Dir filesep 'lifeTime_velocity.fig']);


% normalize speed distribution

% plot speed histogram
[h,x] = hist(track_data(:,4),10);
h_area =  trapz(x,h);
h = h./ h_area; 
figure
plot(x,h);
xlabel('Speckle speed (pixels/frame)');

% normalize life time distribution

% plot life time histogram
clear x;
clear h;
x=2:60;
[h,x] = hist(track_data(:,3),x);
h_area =  trapz(x,h);
h = h./ h_area; 
figure
plot(x,h);
xlabel('Speckle life time');


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

do_segmentation = 1;

% select region to clalculate LP/ LA ratio
%[name, path] = uigetfile('*','Select reference image for segment overlay');
%img = imread([path, name]);
%[BW,x_region,y_region] = roipoly(img);       


while do_segmentation
    output = cytoTrackAnalysisDlg;

    life_time_thresh = output.life_time_thresh;
    track_data_short = track_data(track_data(:,3) <= life_time_thresh,:);
    track_data_short = track_data_short(track_data_short(:,3) >  2,:);
    track_data_long  = track_data(track_data(:,3) > life_time_thresh, :);

    vel_thresh  = output.speed_thresh;
    track_data_slow  = track_data(track_data(:,4) <= vel_thresh,:);
    track_data_fast  = track_data(track_data(:,4) >  vel_thresh, :);
 
    track_data_fast_long = track_data_fast(track_data_fast(:,3) > life_time_thresh, :);

        
    h = figure;
    plot(track_data_fast(:,2),track_data_fast(:,1),'oc','MarkerFaceColor','c','MarkerSize',2);
    %set(); 
    hold on
    plot(track_data_long(:,2),track_data_long(:,1),'ob','MarkerFaceColor','b','MarkerSize',2);
    title(['Fast tracks. Speed > ' num2str(vel_thresh) ' (pixel/frames)' '   Long tracks. Life time > ' num2str(life_time_thresh) ' (frames)']);
    axis off
    axis ij
    axis equal
    hgsave(h,[result_Dir filesep 'lp_la_segmentation.fig']); 
    print(h, [result_Dir filesep 'lp_la_segmentation.eps'],'-depsc2','-tiff'); 
    print(h, [result_Dir filesep 'lp_la_segmentation.tif'],'-dtiff');     
    
    
    % analyze speckle populations
    % total number of speckles
    total_num       = size(track_data,1);
    lp_num          = size(track_data_fast,1);
    la_num          = size(track_data_long,1);
    la_lp_num       = size(track_data_fast_long,1);
    total_num;
    percent_lp          = lp_num/total_num;
    percent_la          = la_num/total_num;
    percent_lp_from_la  = la_lp_num/la_num;
    lp_to_la  = lp_num/la_num;
    
    disp(['LP/LA        = ', num2str(lp_to_la)]);
    disp(['LP/Total     = ', num2str(percent_lp)]); 
    disp(['LA/Total     = ', num2str(percent_la)]);
    disp(['LP in LA     = ', num2str(percent_lp_from_la)]);
        
    % Average speed of fast speckles
    av_speed_fast_speckles = mean(track_data_fast(:,4));
    
    % Average speed of stable speckles
    av_speed_stable_speckles = mean(track_data_long(:,4));
    
    disp(['Average speed fast speckles =   ', num2str(av_speed_fast_speckles)]);
    disp(['Average speed stable speckles = ', num2str(av_speed_stable_speckles)]);
    fprintf(1,'\n\n');
    
    % save all data
    save([result_Dir filesep 'track_data'],'track_data');
    slow_sores = ones(size(track_data_slow));
    fast_sores = ones(size(track_data_fast));
    slow_scores(:,2:3) = track_data_slow(:,1:2);
    slow_scores(:,4)   = track_data_slow(:,4);
    save([result_Dir filesep 'slow_scores'],'slow_scores');
    fast_scores(:,2:3) = track_data_fast(:,1:2);
    fast_scores(:,4)   = track_data_fast(:,4);
    save([result_Dir filesep 'fast_scores'],'fast_scores');


%     t_n
%     [hi x_bin] = hist(track_length,40);
%     figure
%     plot(x_bin,hi);
 
    % dialog to ask if user likes to repeat the segmentation
    button_ans = questdlg('Repeat segmentation?');
    if strcmp(button_ans, 'Yes')
        do_segmentation = 1;
        clear slow_scores;
        clear fast_scores;
    else
        do_segmentation = 0;        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map = zeros(max(track_data(:,1)),max(track_data(:,2)));

for i=1:size(track_data_fast,1)
    map(int16(track_data_fast(i,1)),int16(track_data_fast(i,2)),1) = 1;
    map(int16(track_data_fast(i,1)),int16(track_data_fast(i,2)),2) = 0;
    map(int16(track_data_fast(i,1)),int16(track_data_fast(i,2)),3) = 0;
end
for i=1:size(track_data_long,1)
    map(track_data_long(i,1),track_data_long(i,2),1) = 0.2;
    map(track_data_long(i,1),track_data_long(i,2),2) = 0.5;
    map(track_data_long(i,1),track_data_long(i,2),3) = 0.8;
end


do_region = 1;
while do_region

    [BW,x_region,y_region] = roipoly(map);

    lp_region = track_data_fast(inpolygon(track_data_fast(:,1),track_data_fast(:,2),y_region,x_region),:);
    la_region = track_data_long(inpolygon(track_data_long(:,1),track_data_long(:,2),y_region,x_region),:);

    lp_num_region = size(lp_region,1); 
    la_num_region = size(la_region,1); 
    
    lp_to_la_region = lp_num_region/la_num_region;
    disp(['LP/LA region = ', num2str(lp_to_la_region)]);
    button_ans = questdlg('Repeat?');
    if strcmp(button_ans, 'Yes')
        do_region = 1;
        clear lp_region;
        clear la_region;
    else
        do_region = 0;
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');



%    track_map_slow = zeros(max(track_data(:,1)),max(track_data(:,2)));
%    track_map_fast = zeros(max(track_data(:,1)),max(track_data(:,2)));
% 
%         for i=1:size(track_data_short,1)
%             track_map_short(int16(track_data_short(i,1)),int16(track_data_short(i,2)),1) = 1;
%             track_map_short(int16(track_data_short(i,1)),int16(track_data_short(i,2)),2) = 0;
%             track_map_short(int16(track_data_short(i,1)),int16(track_data_short(i,2)),3) = 0;
%         end
%         for i=1:size(track_data_long,1)
%             track_map_long(track_data_long(i,1),track_data_long(i,2),1) = 1;
%             track_map_long(track_data_long(i,1),track_data_long(i,2),2) = 0;
%             track_map_long(track_data_long(i,1),track_data_long(i,2),3) = 1;
%         end
%         figure,imshow(track_map_short,[]);
%         title(['Short tracks. Life time <= ' num2str(threshold) ' frames']);
%         imwrite(track_map_short,[result_Dir filesep 'track_map_short.tif'],'tif');
%         figure,imshow(track_map_long,[]);
%         title(['Long tracks. Life time > ' num2str(threshold) ' frames']);
%         imwrite(track_map_long,[result_Dir filesep 'track_map_long.tif'],'tif');
% 
%    
%         for i=1:size(track_data_slow,1)
%             track_map_slow(int16(track_data_slow(i,1)),int16(track_data_slow(i,2)),1) = 1;
%             track_map_slow(int16(track_data_slow(i,1)),int16(track_data_slow(i,2)),2) = 0;
%             track_map_slow(int16(track_data_slow(i,1)),int16(track_data_slow(i,2)),3) = 0;
%         end
%         for i=1:size(track_data_fast,1)
%             track_map_fast(track_data_fast(i,1),track_data_fast(i,2),1) = 1;
%             track_map_fast(track_data_fast(i,1),track_data_fast(i,2),2) = 0;
%             track_map_fast(track_data_fast(i,1),track_data_fast(i,2),3) = 1;
%         end
%         figure,imshow(track_map_slow,[]);
%         title(['Slow tracks. Speed <= ' num2str(vel_thresh) ' (pixel/frames)']);
%         imwrite(track_map_slow,[result_Dir filesep 'track_map_slow.tif'],'tif');
%         figure,imshow(track_map_fast,[]);
%         title(['Fast tracks. Speed > ' num2str(vel_thresh) ' (pixel/frames)']);
%         imwrite(track_map_fast,[result_Dir filesep 'track_map_fast.tif'],'tif');

