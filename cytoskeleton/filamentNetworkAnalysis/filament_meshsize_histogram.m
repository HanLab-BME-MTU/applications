function [network_feature_channel_frame, network_feature_channel] = filament_meshsize_histogram(MD, ROI, radius,pace,range)
% function to do network analysis with input MD

% input:    MD:    the loaded movieData object.
%           ROI:   the user input ROI, if not defined [], will use the whole
%                  image area
%           radius: in the case of no cell segmentation is given, build a
%                   cell area by imdilate of this radius
% output:   network_feature, a cell structure for each channel, each frame.
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

%
%% find packages and processes indices
movie_Dir = MD.outputDirectory_;

% find all the index of different processes
package_process_ind_script;
network_feature_channel_frame=cell(length(MD.channels_),nFrame);
network_feature_channel=cell(length(MD.channels_),1);

outdir = [MD.outputDirectory_,filesep,'meshsize_results'];

% make out put directory if not existing
if(~exist(outdir,'dir'))
    mkdir(outdir);
end

%% Analysis
% for each channel do analysis
for iChannel = 1 :  length(MD.channels_)
    channel_meshsize_local_maximum_pool=[];
    
    % for each frame
    for iFrame = 1 : nFrame
        display(['iChannel: ', num2str(iChannel),', iFrame:', num2str(iFrame)]);
        meshsize_local_maximum_pool=[];
        %% % Load the data
        % load the filament segmentation results
        VIF_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
        VIF_current_seg = isnan(VIF_orientation)==0;
        
        % if the input ROI is [], then use the whole area
        if(isempty(ROI))
            ROI = ones(size(VIF_current_seg));
            
            if(indexCellRefineProcess>0)
                Cell_Mask = ROI.*((MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame))>0);
            else
                Cell_Mask = ROI.*(imdilate(VIF_current_seg,ones(2*radius+1,2*radius+1)));
            end
        else
            % if there is input ROI, we still want to exclude the part not
            % in cell if segmentation is given.
            
            if(indexCellRefineProcess>0)
                Cell_Mask = ROI.*(MD.processes_{indexCellRefineProcess}.loadChannelOutput(iChannel,iFrame));
            else
                % if ROI is given, but no cell segmentation, just use ROI
                Cell_Mask = ROI;
            end
        end
        output_feature = filament_bw_meshsize_histogram(VIF_current_seg, Cell_Mask, radius,pace,range);
        
        network_feature_channel_frame{iChannel,iFrame} = output_feature;
        
        % stack for the all the frames in this movie
        channel_meshsize_local_maximum_pool = [channel_meshsize_local_maximum_pool;  output_feature.meshsize_local_maximum_pool];
    end
    
    
    [h,bin]= hist(channel_meshsize_local_maximum_pool,0:pace:range);
    h = h./(sum(h))*100;
    
    % find the mode
    ind = find(h==max(h));
    ind = ind(1);
    mode_bin = bin(ind);
    mode_h = h(ind);
    
    % find the mean
    mean_bin = double(mean(channel_meshsize_local_maximum_pool));
    
    h4 = figure(4); hold off;
    
    bar(bin, h);
    
    axis([0 range+pace 0 max(h)+10]);
    
    real_axis=  axis;
    hold on;
    plot(mode_bin, mode_h, 'r*');
    text(mode_bin-2, mode_h+1, ['Mode: ',num2str(mode_bin)]);
    
    plot([mean_bin mean_bin],[0 real_axis(4)],'m');
    
    text(mean_bin, real_axis(4)-1.5, ['Mean: ',num2str(mean_bin)]);
    
    title(['Meshsize Measurement, Channel:',num2str(iChannel),', all Frames']);
    xlabel('distance to filament (unit: pixel)');
    ylabel('Percentage(%)');
    
    saveas(h4,[outdir, filesep, 'dist_hist_ch',num2str(iChannel),'_allframe.tif']);
    
    % save the results
    output_feature = [];
    output_feature.channel_meshsize_local_maximum_pool = channel_meshsize_local_maximum_pool;
    output_feature.mode_meshsize = mode_bin;
    output_feature.mean_meshsize = mean_bin;
    
    network_feature_channel{iChannel,1} = output_feature;
    
    
end


save([outdir, filesep, 'dist_hist_output.mat'],'network_feature_channel_frame','network_feature_channel');
