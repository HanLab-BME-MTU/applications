function filament_network_plotting_with_2MD(MD1, MD2, varargin)
% function to plot network with the convenience of control parameters
% Input:
%   MD1, MD2: the two loaded movieData objects
% Optional Input:
%   preset_index:  the fastest way to get plots, with index of 1, both 
%                   channels no 1, red_green_blackbackground; or index of 2,
%                   both channels no 1, red_green_whitebackground. Default, 1
%   preset:        the string way to lookup for preset,
%
%   'index_channels': the channels to be plotted in two MD1s, default [1 1]
%   'network_color_cells': cell for the color vectors for each channels, in the form of three values for RGB
%   'background_color': the background color vector
%   'frames_index': array for frame index, default, all frames
%   'line_widths': array for the linewidth, if line, marker size if marker;
%   'overlay': flag as to if overlay the filament network on image or background
%   'image_inversion': flag, in the case of single channel image display, flag to say if invert the color or not
%   'int_max': the maximum in inversion of greylevel image, e.g., 255
%   'background_type':  'solid_color' or 'image'
%   'orientation_heatmap': flag, if orientation heatmap is asked for, put 1 in the vector
%   'marker_type': the marker type for filament display, e.g.,'.',':'...
%   'pixel_based': flag for pixel based method, if 1, construct pixel RGB images
%   'pixel_dilate_size': for pixel based method, dilate the pixel with this size
%   'image_background_ch': the grey image channel as background
%   'axis_limit': a 4 number vector to limit the showing window, like [ 1 200 1 400]

%  Test try examples:
%
%  filament_network_plotting_with_nMD1(MD11, MD12)
%           to get default plotting of two channel 1 white background and red
%           and green lines 

%  filament_network_plotting_with_nMD1(MD11, MD12,'pixel_based',1,'pixel_dilate_size',2,'image_background_ch',1,'image_inversion',1)
%           to plot pixel based version, with dilation of 2, with
%           background of image channel 1 
%
%  filament_network_plotting_with_nMD1(MD11, MD12,'marker_type',{'d','.'},'line_widths',[20 10])
%           to plot with default plotting with different marker and size
%
%  filament_network_plotting_with_nMD1(MD11, MD12,'pixel_based',1,'pixel_dilate_size',3,'overlay,1,'background_type','image')
%           to plot the original RGB tiff image, and overlay with filaments
%
%  filament_network_plotting_with_nMD1(MD11, MD12,'index_channels',[2 1],'network_color_cells',{[1 0 0],[ 0 1 0]},'overlay',1,'background_type','image','image_inversion',1,'image_background_ch',1)
%           to plot the channel 1 image in grey and overlay with filament network, with user defined colors
%
%
ip = inputParser;
ip.addRequired('MD1',@(x)(isa(x,'MovieData')));
ip.addOptional('preset_index', 1, @isnumeric)
preset_pool = {'red_green_blackbackground', 'red_green_whitebackground'};
ip.addOptional('preset_name', 'red_green_whitebackground', @(x) ismember(x, preset_pool));

ip.addOptional('index_channels', 0, @isvector);
ip.addOptional('network_color_cells', {[ 0 1 0], [1 0 0]}, @iscell);
ip.addOptional('background_color', [0 0 0], @isvector);
ip.addOptional('frames_index', [1 2], @isvector);
ip.addOptional('line_widths', [2 2], @isvector);
ip.addOptional('overlay', 1, @isscalar);
ip.addOptional('image_inversion', 0, @isscalar);
ip.addOptional('image_background_ch', 1, @isscalar);
ip.addOptional('int_max', 255, @isscalar);

background_pool = {'image', 'solid_color'};
ip.addOptional('background_type', 'solid_color', @(x) ismember(x, background_pool));
ip.addOptional('orientation_heatmap', [0 0], @isvector);
ip.addOptional('marker_type', {'.','.'}, @iscell);
ip.addOptional('pixel_based', 1, @isscalar);
ip.addOptional('pixel_dilate_size', 1, @isscalar);
ip.addOptional('axis_limit', [1 1 1 1], @vector);

all_options = {'preset_index','preset_name','index_channels',...
    'network_color_cells','background_color','frames_index',...
    'line_widths','overlay','image_inversion','int_max',...
    'background_type','orientation_heatmap','marker_type',...
    'pixel_based','pixel_dilate_size','image_background_ch',...
    'axis_limit'};

ip.parse(MD1, varargin{:});

input_setting =ip.Results;

% index_channels
% frames_index
% background_color
%
% network_color_cells

nFrames = MD1.nFrames_;

% Set all the default parameter set, just as in preset #1 and "red_green_whitebackground"
param =[];
param.index_channels = [1 1];
param.background_color = [1 1 1];
param.network_color_cells =  {[0 1 0], [1 0 0]};
param.frames_index = 1:nFrames;
param.line_widths =  [2 2];
param.overlay =  1;
param.image_inversion = 0;
param.int_max =  255;
param.background_type = 'solid_color';
param.orientation_heatmap = [0 0];
param.marker_type = {'.','.'};
param.pixel_based = 0;
param.pixel_dilate_size = 1;
param.image_background_ch = 0;
param.axis_limit = [1 MD1.imSize_(2) 1  MD1.imSize_(1)];

user_assigned_flag = ones(1,20);


for iTerms = 1 : length(all_options)
    if(ismember(all_options{iTerms},ip.UsingDefaults));
        user_assigned_flag(iTerms) = 0;
    end
end

% if preset_index is 2, then black background
if (user_assigned_flag(1) == 1 && input_setting.preset_index == 2)
    param.background_color = [0 0 0];
end

% if preset_name ask for black background, then black background
% as you can see, if preset_index ==1, and preset_name ask for black,
% then preset_name wins
if (user_assigned_flag(3) == 1)
    param.index_channels = input_setting.index_channels;
end

% again, if assigned directly by user, this over write the previous preset
if (user_assigned_flag(4) == 1)
    param.network_color_cells = input_setting.network_color_cells;
end

if (user_assigned_flag(5) == 1)
    param.background_color = input_setting.background_color;
end

if (user_assigned_flag(6) == 1)
    param.frames_index = input_setting.frames_index;
end


if (user_assigned_flag(7) == 1)
    param.line_widths = input_setting.line_widths;
end

if (user_assigned_flag(8) == 1)
    param.overlay = input_setting.overlay;
end

if (user_assigned_flag(9) == 1)
    param.image_inversion = input_setting.image_inversion;
end

if (user_assigned_flag(10) == 1)
    param.int_max = input_setting.int_max;
end

if (user_assigned_flag(11) == 1)
    param.background_type = input_setting.background_type;
end

if (user_assigned_flag(12) == 1)
    param.orientation_heatmap = input_setting.orientation_heatmap;
end

if (user_assigned_flag(13) == 1)
    param.marker_type = input_setting.marker_type;
end

if (user_assigned_flag(14) == 1)
    param.pixel_based = input_setting.pixel_based;
end

if (user_assigned_flag(15) == 1)
    param.pixel_dilate_size = input_setting.pixel_dilate_size;
end

if (user_assigned_flag(16) == 1)
    param.image_background_ch = input_setting.image_background_ch;
end

if (user_assigned_flag(17) == 1)
    param.axis_limit = input_setting.axis_limit;
end

% sanity checking
nChannels =  length(param.index_channels);
if(length(param.network_color_cells)~=nChannels)
    msgbox('Network colors cell defined is invalid');
    return;
end

if( param.frames_index> nFrames)
    msgbox('Frame index defined is invalid');
    return;
end

if(length(param.line_widths)~=nChannels)
    msgbox('Line widths vector defined is invalid');
    return;
end

if(param.image_inversion==1 && param.image_background_ch==0)
    msgbox('To show inverted image,define the channel');
    return;
end

if(param.int_max<0)
    msgbox('Please enter valid maximum intensity.');
    return;
end

if(param.axis_limit(1)<0 || param.axis_limit(2)> MD1.imSize_(2) ...
        || param.axis_limit(3)<0 || param.axis_limit(4)> MD1.imSize_(1))
    msgbox('Please enter valid axis limit.');
    return;
end

% start plotting

movie_Dir = MD1.outputDirectory_;

outdir = [movie_Dir,filesep,'network_plots'];
if(~exist(outdir,'dir'))
    mkdir(outdir);
end

package_process_ind_script;

for iFrame = param.frames_index
    display(['Plotting for frame: ', num2str(iFrame)]);
    
    % if pixel based plot is asked for
    if(param.pixel_based>0)
        
        RGB_channels_seg = zeros(MD1.imSize_(1),MD1.imSize_(2),3);
        
        RGB_channels_img = zeros(MD1.imSize_(1),MD1.imSize_(2),3);
        RGB_channels_img(:,:,1) = param.background_color(1);
        RGB_channels_img(:,:,2) = param.background_color(2);
        RGB_channels_img(:,:,3) = param.background_color(3);
        
        RGB_channels_user_color = zeros(MD1.imSize_(1),MD1.imSize_(2),3);
        R_channels_user_color = RGB_channels_user_color(:,:,1);
        G_channels_user_color = RGB_channels_user_color(:,:,2);
        B_channels_user_color = RGB_channels_user_color(:,:,3);
        
        for iChInd = 1: nChannels
            iChannel = param.index_channels(iChInd);
            
            VIF_orientation = MD1.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
            VIF_current_model = MD1.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
            
            VIF_current_seg = (isnan(VIF_orientation)==0);
            
            VIF_current_seg = iMD1ilate(VIF_current_seg, ...
                ones(param.pixel_dilate_size,param.pixel_dilate_size));
            
            % if flattened use flatten, otherwise use original
            if indexFlattenProcess~=0
                VIF_img =  MD1.processes_{indexFlattenProcess}.loadChannelOutput(iChannel,iFrame);
            else
                VIF_img =  MD1.channels_(iChannel).loadImage(iFrame);
            end
            
            if(param.image_inversion==1)
                VIF_img = param.int_max-VIF_img;
            end
            
            RGB_channels_img(:,:,iChInd) = VIF_img;
            
            RGB_channels_seg(:,:,iChInd) = VIF_current_seg;
            
            R_channels_user_color = R_channels_user_color + VIF_current_seg*param.network_color_cells{iChInd}(1);
            G_channels_user_color = G_channels_user_color + VIF_current_seg*param.network_color_cells{iChInd}(2);
            B_channels_user_color = B_channels_user_color + VIF_current_seg*param.network_color_cells{iChInd}(3);
            
            if(param.image_inversion==1 && param.image_background_ch==iChannel)
                RGB_channels_img(:,:,1) = VIF_img;
                RGB_channels_img(:,:,2) = VIF_img;
                RGB_channels_img(:,:,3) = VIF_img;
            end
            
            if(param.image_background_ch==iChannel)
                RGB_channels_img(:,:,1) = VIF_img;
                RGB_channels_img(:,:,2) = VIF_img;
                RGB_channels_img(:,:,3) = VIF_img;
            end
            
            
        end
        
        
        h5=figure(5);imagesc(RGB_channels_img/255);axis equal;axis off;
        axis(param.axis_limit);
        set(gcf,'Units','normal');
        set(gca,'Position',[0 0 1 1]);
        iptsetpref('ImshowBorder','tight');
        saveas(h5,[outdir,filesep,'RGB_img_frame_',num2str(iFrame),'.tif']);
        saveas(h5,[outdir,filesep,'RGB_img_frame_',num2str(iFrame),'.fig']);
        print(h5,'-depsc',[outdir,filesep,'RGB_img_frame_',num2str(iFrame),'.eps']);
        
        % change back ground color
        all_channel_current_seg = sum(RGB_channels_seg,3)>0;
        
        % these are for the network alone
        if(strcmp(param.background_type,'solid_color'))
            
            % if there is no segmentation in any channel, set it to background color
            R_channels_user_color(find(all_channel_current_seg==0)) = param.background_color(1);
            G_channels_user_color(find(all_channel_current_seg==0)) = param.background_color(2);
            B_channels_user_color(find(all_channel_current_seg==0)) = param.background_color(3);
            
            R_channels_user_color(R_channels_user_color>1)=1;
            G_channels_user_color(G_channels_user_color>1)=1;
            B_channels_user_color(B_channels_user_color>1)=1;
            
            RGB_channels_user_color(:,:,1) = double(R_channels_user_color);
            RGB_channels_user_color(:,:,2) = double(G_channels_user_color);
            RGB_channels_user_color(:,:,3) = double(B_channels_user_color);
        else
            % if there is no segmentation in any channel, set it to single image
            if( param.overlay>0)
                R_channels_img = RGB_channels_img(:,:,1);
                G_channels_img = RGB_channels_img(:,:,2);
                B_channels_img = RGB_channels_img(:,:,3);
                
                R_channels_user_color(find(all_channel_current_seg==0)) = R_channels_img(find(all_channel_current_seg==0))/255;
                G_channels_user_color(find(all_channel_current_seg==0)) = G_channels_img(find(all_channel_current_seg==0))/255;
                B_channels_user_color(find(all_channel_current_seg==0)) = B_channels_img(find(all_channel_current_seg==0))/255;
                
                R_channels_user_color(R_channels_user_color>1)=1;
                G_channels_user_color(G_channels_user_color>1)=1;
                B_channels_user_color(B_channels_user_color>1)=1;
                
                RGB_channels_user_color(:,:,1) = double(R_channels_user_color);
                RGB_channels_user_color(:,:,2) = double(G_channels_user_color);
                RGB_channels_user_color(:,:,3) = double(B_channels_user_color);
            end
        end
        
        if(strcmp(param.background_type,'solid_color') || ...
                (~strcmp(param.background_type,'solid_color') && param.overlay) )
            h6 = figure(6);imagesc(RGB_channels_user_color(:,1:end,:));axis equal;axis off;
            figure(6);
            axis(param.axis_limit);
            %             set(gca, 'Position', get(gca, 'OuterPosition') - ...
            %                 get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
            set(gcf,'Units','normal');
            set(gca,'Position',[0 0 1 1]);
            iptsetpref('ImshowBorder','tight');
            saveas(h6,[outdir,filesep,'RGB_seg_frame_',num2str(iFrame),'.fig']);
            saveas(h6,[outdir,filesep,'RGB_seg_frame_',num2str(iFrame),'.tif']);
            print(h6,'-depsc',[outdir,filesep,'RGB_seg_frame_',num2str(iFrame),'.eps']);
        end
        
    else
        
        % figure 3 is for direct plot with the filament network model obtained
        background_RGB = zeros(MD1.imSize_(1),MD1.imSize_(2),3);
        
        % first if solid color back ground is choosen
        if(strcmp(param.background_type,'solid_color'))
            background_RGB = zeros(MD1.imSize_(1),MD1.imSize_(2),3);
            background_RGB(:,:,1) = param.background_color(1)*255;
            background_RGB(:,:,2) = param.background_color(2)*255;
            background_RGB(:,:,3) = param.background_color(3)*255;
        else
            % second, if image back ground is choosen
            % Get the background ready without loading the network
            for iChInd = 1: nChannels
                iChannel = param.index_channels(iChInd);
                
                % if flattened use flatten, otherwise use original
                if indexFlattenProcess~=0
                    VIF_img =  MD1.processes_{indexFlattenProcess}.loadChannelOutput(iChannel,iFrame);
                else
                    VIF_img =  MD1.channels_{iChannel}.loadImage(iFrame);
                end
                
                % if user ask for inversion
                if(param.image_inversion==1)
                    VIF_img = param.int_max - VIF_img;
                    % if the user just want to show this channel
                    if(param.image_background_ch==iChannel && strcmp(param.background_type,'image'))
                        background_RGB(:,:,1) = VIF_img;
                        background_RGB(:,:,2) = VIF_img;
                        background_RGB(:,:,3) = VIF_img;
                        break;
                    end
                end
                
                % if user didn't ask for inversion, but just want this channel in
                % gray color
                if(param.image_background_ch==iChannel && strcmp(param.background_type,'image'))
                    background_RGB(:,:,1) = VIF_img;
                    background_RGB(:,:,2) = VIF_img;
                    background_RGB(:,:,3) = VIF_img;
                    break;
                end
                
                % if the user want all channels to be displayed in color
                if(param.image_background_ch==0 && strcmp(param.background_type,'image'))
                    background_RGB(:,:,iChInd) = VIF_img;
                end
            end
        end
        
        h4 = figure(4); hold off;
        imagesc(background_RGB(:,1:end,:)/255);axis equal;axis off;
        hold on;
        
        % if background is solid color, don't save anything, else, save image.
        if(~strcmp(param.background_type,'solid_color'))
            axis(param.axis_limit);
            %             set(gca, 'Position', get(gca, 'OuterPosition') - ...
            %                 get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
            set(gcf,'Units','normal');
            set(gca,'Position',[0 0 1 1]);
            
            saveas(h4,[outdir,filesep,'image_frame_',num2str(iFrame),'.tif']);
            saveas(h4,[outdir,filesep,'image_frame_',num2str(iFrame),'.fig']);
            print(h4,'-depsc',[outdir,filesep,'image_frame_',num2str(iFrame),'.eps']);
        end
        
        % if segmentation overlay is requested
        if(param.overlay==1)
            for iChInd = 1: nChannels
                iChannel = param.index_channels(iChInd);
                
                VIF_orientation = MD1.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
                VIF_current_model = MD1.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
                VIF_heat_output = MD1.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','RGB_seg_orient_heat_map');
                
                
                [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
                    = filament_model_to_digital_with_orientation(VIF_current_model);
                
                
                if(param.orientation_heatmap(iChInd)==1)
                    h_iChInd = figure(iChInd);
                    imagesc(VIF_heat_output);
                    axis equal;axis off;
                    axis(param.axis_limit);
                    
                    %                     set(gca, 'Position', get(gca, 'OuterPosition') - ...
                    %                         get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                    set(gcf,'Units','normal');
                    set(gca,'Position',[0 0 1 1]);
                    
                    saveas(h_iChInd,[outdir,filesep,'image_frame_',num2str(iFrame),'.tif']);
                    saveas(h_iChInd,[outdir,filesep,'image_frame_',num2str(iFrame),'.fig']);
                    print(h_iChInd,'-depsc',[outdir,filesep,'image_frame_',num2str(iFrame),'.eps']);
                end
                
                h4 = figure(4);
                if(strcmp(param.marker_type{iChInd},'-') || ...
                        strcmp(param.marker_type{iChInd},'--') || ...
                        strcmp(param.marker_type{iChInd},':') || ...
                        strcmp(param.marker_type{iChInd},'-.') )
                    for iM = 1 : length(Vif_digital_model)
                        hold on;
                        plot(Vif_digital_model{iM}(:,1),Vif_digital_model{iM}(:,2),...
                            '-',...
                            'Color',param.network_color_cells{iChInd}',...
                            'LineWidth',param.line_widths(iChInd));
                    end
                else
                    plot(VIF_XX,VIF_YY,...
                        param.marker_type{iChInd},...
                        'Color',param.network_color_cells{iChInd}',...
                        'MarkerSize',param.line_widths(iChInd));
                end
            end
            
            h4=figure(4);
            axis(param.axis_limit);
            %             set(gca, 'Position', get(gca, 'OuterPosition') - ...
            %                 get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
            %
            
            set(gcf,'Units','normal');
            set(gca,'Position',[0 0 1 1]);
            iptsetpref('ImshowBorder','tight');
            
            % set(gcf, 'PaperUnits','centimeter');
            % set(gcf, 'PaperSize', [MD1.imSize_(2)/20 MD1.imSize_(1)/20]);
            % set(gcf, 'PaperPositionMode', 'manual');
            % set(gcf, 'PaperPosition',[0 0  MD1.imSize_(2)/20 MD1.imSize_(1)/20]);
            % set(gcf, 'Position',[100 50  MD1.imSize_(2) MD1.imSize_(1)]);
            %
            % set(gca,'DataAspectRatioMode','auto');
            % set(gca,'PlotBoxAspectRatioMode','auto');
            % set(gca,'CameraViewAngleMode','auto');
            %
            % set(gca,'Units','pixels');
            % set(gca, 'Position',[0 0  MD1.imSize_(2) MD1.imSize_(1)]);
            
            saveas(h4,[outdir,filesep,'filament_overlay_',num2str(iFrame),'.tif']);
            saveas(h4,[outdir,filesep,'filament_overlay_',num2str(iFrame),'.fig']);
            print(h4,'-depsc',[outdir,filesep,'filament_overlay_',num2str(iFrame),'.eps']);
            
            
            
            
        end
    end
end

