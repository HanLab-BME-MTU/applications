function network_plotting_with_MD(MD, varargin)
% function to plot network with the convenience of control parameters
% Input:
%   MD: the loaded movieData object
% Optional Input:
%   preset_index:  the fastest way to get plots, with index of 1, two
%                   channels, red_green_blackbackground; or index of 2,
%                   two channels, red_green_whitebackground. Default, 1
%   preset:        the string way to lookup for preset,
%
%   index_channels: the channels to be plotted, default [1 2]
%   network_color_cells: color vectors for each channels, in the form of three values for RGB
%   background_color: the background color vector
%   frames_index: array for frame index, default, all frames

ip =inputParser;
ip.addRequired('MD',@(x)(isa(x,'MovieData')));
ip.addOptional('preset_index', 1, @isnumeric)
preset_pool = {'red_green_blackbackground', 'red_green_whitebackground'};
ip.addOptional('preset_name', 'red_green_whitebackground', @(x) ismember(x, preset_pool));

ip.addOptional('index_channels', 0, @isvector);
ip.addOptional('network_color_cells', {[ 0 1 0], [1 0 0]}, @iscell);
ip.addOptional('background_color', [0 0 0], @isvector);
ip.addOptional('frames_index', [1 2], @isvector);
ip.addOptional('line_widths', [2 2], @isvector);
ip.addOptional('overlay', 0, @isscalar);
ip.addOptional('image_inversion', 1, @isscalar);

all_options = {'preset_index','preset_name','index_channels',...
    'network_color_cells','background_color','frames_index',...
    'line_widths','overlay','image_inversion'};

ip.parse(MD, varargin{:});

input_setting =ip.Results;

% index_channels
% frames_index
% background_color
%
% network_color_cells

nFrames = MD.nFrames_;

% Set all the default parameter set, just as in preset #1 and "red_green_whitebackground"
param =[];
param.index_channels = [1 2];
param.background_color = [1 1 1];
param.network_color_cells =  {[ 1 0 0], [0 1 0]};
param.frames_index = 1:nFrames;
param.line_widths =  [2 2];
param.overlay =  0;
param.image_inversion = 1;

user_assigned_flag = ones(1,7);

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

% start plotting

movie_Dir = MD.outputDirectory_;

outdir = [movie_Dir,filesep,'network_plots'];
mkdir(outdir);

package_process_ind_script;

for iFrame = 1 : nFrames
    display(['Plotting for frame: ', num2str(iFrame)]);
    RGB_channels_seg = zeros(MD.imSize_(1),MD.imSize_(2),3);
    RGB_channels_img = zeros(MD.imSize_(1),MD.imSize_(2),3);
    RGB_channels_user_color = zeros(MD.imSize_(1),MD.imSize_(2),3);
    R_channels_user_color = zeros(MD.imSize_(1),MD.imSize_(2));
    G_channels_user_color = zeros(MD.imSize_(1),MD.imSize_(2));
    B_channels_user_color = zeros(MD.imSize_(1),MD.imSize_(2));
    
    if(param.image_inversion==1)
        RGB_channels_img(:)=255;
    end
    
    h3=figure(3);
    
    % figure 3 is for direct plot with the filament network model obtained
   
    
    
    for iChInd = 1: nChannels
        iChannel = param.index_channels(iChInd);
        
        VIF_orientation = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
        VIF_current_model = MD.processes_{indexFilamentSegmentationProcess}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
        
        [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
            = filament_model_to_digital_with_orientation(VIF_current_model);
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        
        VIF_img =  MD.processes_{3}.loadChannelOutput(1,iFrame);
        
        if(param.image_inversion==1)
            VIF_img = 255-VIF_img;
        end
        
        RGB_channels_img(:,:,iChInd) = VIF_img;
        RGB_channels_seg(:,:,iChInd) = VIF_current_seg;
        
        R_channels_user_color = R_channels_user_color + VIF_current_seg*param.network_color_cells{iChInd}(1);
        G_channels_user_color = G_channels_user_color + VIF_current_seg*param.network_color_cells{iChInd}(2);
        B_channels_user_color = B_channels_user_color + VIF_current_seg*param.network_color_cells{iChInd}(3);
        
    end
    
    % if image overlay is requested
    if(param.overlay==1)
        h1=figure(1);imagesc(RGB_channels_img/255);axis equal;axis off;
        saveas(h1,[outdir,filesep,'RGB_img_frame_',num2str(iFrame),'.tif']);
        saveas(h1,[outdir,filesep,'RGB_img_frame_',num2str(iFrame),'.fig']);
    end
    
    % these are for the network alone
    
    % change back ground color
    all_channel_current_seg = sum(RGB_channels_seg,3)>0;
    % if there is no segmentation in any channel, set it to background color
    R_channels_user_color(find(all_channel_current_seg==0)) = param.background_color(1);
    G_channels_user_color(find(all_channel_current_seg==0)) = param.background_color(2);
    B_channels_user_color(find(all_channel_current_seg==0)) = param.background_color(3);
    
    RGB_channels_user_color(:,:,1) = double(R_channels_user_color);
    RGB_channels_user_color(:,:,2) = double(G_channels_user_color);
    RGB_channels_user_color(:,:,3) = double(B_channels_user_color);
    
    h2 = figure(2);imagesc(RGB_channels_user_color(:,1:end,:));axis equal;axis off;
    figure(2);
    %     set(gca, 'Position', get(gca, 'OuterPosition') - ...
    %       get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    saveas(h2,[outdir,filesep,'RGB_seg_frame_',num2str(iFrame),'.tif']);
    saveas(h2,[outdir,filesep,'RGB_seg_frame_',num2str(iFrame),'.fig']);
   
    
     
end

