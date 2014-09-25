function output_feature = filament_bw_meshsize_histogram(VIF_current_seg, Cell_Mask, radius,pace,range)
% function to do network analysis with input MD

% input:    VIF_current_seg:    network segmentation
%           Cell_Mask:   the user input ROI
%           radius: in the case of no cell segmentation is given, build a
%                   cell area by imdilate of this radius
% output:   output_feature, a cell structure for this channel, this frame.
%           Each struct with field of:
%           'straightness_per_filament_pool','orientation_pixel_pool_display',...
%           'length_per_filament_pool'.

%% check radius first
% radius is for if there is no definition of cell area(from segmentation)
% nor is there ROI

% if there is no input requirement for radius or just [], set it as default
if(nargin==2)
    radius = 15;
end
if(nargin==3)
    if(isempty(radius))
        radius = 15;
    end
end


if(nargin<=3)
    pace=3;
end

if(nargin<=4)
    range=36;
end

%% % do mesh size statistics
%         tic
VIF_current_seg = imresize(VIF_current_seg, 1);
Cell_Mask = imresize(Cell_Mask, 1)>0;

dist_map = bwdist(VIF_current_seg);

%         toc

meshsize_map = 2*dist_map;


[cent, varargout]=FastPeakFind(meshsize_map.*imdilate(Cell_Mask, ones(5,5)),2);


maximum_x= cent(1:2:end);
maximum_y= cent(2:2:end);

ind_mask = find(Cell_Mask(sub2ind(size(meshsize_map),maximum_y,maximum_x))>0);

maximum_x = maximum_x(ind_mask);
maximum_y = maximum_y(ind_mask);

h3=figure(3);
imagesc(meshsize_map.*Cell_Mask);
axis image; axis off;
hold on; plot(maximum_x,maximum_y,'m.');

title(['Distmap, Channel:',num2str(iChannel),', Frame:',num2str(iFrame)]);

meshsize_local_maximum_pool = meshsize_map(sub2ind(size(meshsize_map),cent(2:2:end),cent(1:2:end)));

[h,bin]= hist(meshsize_local_maximum_pool,0:pace:range);
h = h./(sum(h))*100;

% find the mode
ind = find(h==max(h));
ind = ind(1);
mode_bin = bin(ind);
mode_h = h(ind);

% find the mean
mean_bin = double(mean(meshsize_local_maximum_pool));
% plot these
h1 =  figure(1); hold off;

bar(bin, h);

axis([0 range+pace 0 max(h)+10]);

real_axis=  axis;
hold on;
plot(mode_bin, mode_h, 'r*');
text(mode_bin-2, mode_h+1, ['Mode: ',num2str(mode_bin)]);

plot([mean_bin mean_bin],[0 real_axis(4)],'m');

text(mean_bin, real_axis(4)-1.5, ['Mean: ',num2str(mean_bin)]);

title(['Meshsize Measurement, Channel:',num2str(iChannel),', Frame:',num2str(iFrame)]);
xlabel('distance to filament (unit: pixel)');
ylabel('Percentage(%)');

saveas(h1,[outdir, filesep, 'dist_hist_ch',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);

% save the results
output_feature = [];
output_feature.meshsize_local_maximum_pool = meshsize_local_maximum_pool;
output_feature.meshsize_map = meshsize_map;
output_feature.meshsize_map = meshsize_map;
output_feature.mode_meshsize = mode_bin;
output_feature.mean_meshsize = mean_bin;
