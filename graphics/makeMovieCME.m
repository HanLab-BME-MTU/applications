% makeMovieCME(data, varargin)
%
% Inputs: 
%             data : data structure from 'loadConditionData()'
%      {trackinfo} : track structure returned by loadtracks()
%
% Options: 
%       'Scalebar' : length of the scalebar. Default: 5 µm
%      'FrameRate' : frame rate of the movie. Default: 15 fps
%           'Zoom' : zoom with respect to pixel resolution of raw data. Default: 1
%

% Francois Aguet (last modified 02/12/2012)

function makeMovieCME(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('tracks', [], @(x) isstruct(x) || isempty(x));
ip.addParamValue('ScaleBar', 5e-6, @isscalar);
ip.addParamValue('FrameRate', 15, @isscalar);
ip.addParamValue('Zoom', 1, @isscalar);
ip.addParamValue('Mode', 'raw', @(x) strcmpi(x, 'raw') | strcmpi(x, 'rgb') | strcmpi(x, 'mask'));
ip.addParamValue('Channel', 1, @(x) ismember(x, 1:length(data.channels)));
ip.addParamValue('FrameRange', 1:data.movieLength);
ip.addParamValue('FileType', 'png', @(x) strcmpi(x, 'png') | strcmpi(x, 'tif') | strcmpi(x, 'tiff'));
ip.addParamValue('DisplayType', 'lifetime', @(x) any(strcmpi(x, {'lifetime', 'category', 'random', 'projection'})));
ip.addParamValue('ShowEvents', false, @islogical);
ip.addParamValue('Detection', []);
ip.addParamValue('ShowGaps', true, @islogical);
ip.addParamValue('FileName', 'Movie', @ischar);
ip.addParamValue('Colormap', []);
ip.parse(data, varargin{:});

% if isunix && ~ismac
%     opengl software; % fix for linux bug (orchestra only?) that causes segmentation fault in 'print'
% end

nx = data.imagesize(2);
ny = data.imagesize(1);
nCh = length(data.channels);
tracks = ip.Results.tracks;
zoom = ip.Results.Zoom;
ch = ip.Results.Channel;
ext = ['.' ip.Results.FileType];
if strcmpi(ip.Results.Mode, 'RGB')
    ch = 1:nCh;
end

mpath = [data.source 'Movies' filesep];
fpath = [mpath 'Frames' filesep];
[~,~] = mkdir(mpath);
[~,~] = mkdir(fpath);


% Compute dynamic range
detectionFile = [data.source 'Detection' filesep 'detection_v2.mat'];
load(detectionFile);
dRange = cell(nCh);
if isfield(frameInfo, 'dRange')
    for c = 1:nCh
        M = arrayfun(@(x) x.dRange{c}, frameInfo, 'UniformOutput', false);
        M = vertcat(M{:});
        dRange{c} = [min(M(1,:)) max(M(2,:))];
    end
else
    for c = 1:nCh
        firstFrame = double(imread(data.framePaths{c}{1}));
        lastFrame = double(imread(data.framePaths{c}{data.movieLength}));
        dRange{c} = [min(min(firstFrame(:)),min(lastFrame(:))) max(max(firstFrame(:)),max(lastFrame(:)))];
    end
end


% Configure figure
h = figure('Visible', 'off', 'Position', [50 50 nx ny]);
iptsetpref('ImshowBorder','tight');
set(h, 'InvertHardcopy', 'off');
set(h, 'PaperUnits', 'Points');
set(h, 'PaperSize', [nx ny]);
set(h, 'PaperPosition', [0 0 nx ny]); % very important
set(h, 'PaperPositionMode', 'auto');
set(h,'DefaultLineLineSmoothing','on'); % points are not rendered !!
set(h,'DefaultPatchLineSmoothing','on');

% Configure axes
ha = axes('Position', [0 0 1 1]);

% Generate frames
fmt = ['%0' num2str(ceil(log10(data.movieLength))) 'd'];
nf = numel(ip.Results.FrameRange);
fprintf('Generating movie frames:     ');
for f = ip.Results.FrameRange
    
    if ~isempty(ip.Results.Detection)
        detection = ip.Results.Detection(f);
    else
        detection = [];
    end
    idx = [tracks.start]<=f & f<=[tracks.end];
    
    if ~isempty(ip.Results.Colormap)
        cmap = ip.Results.Colormap(idx,:);
    else
        cmap = [];
    end
    
    plotFrame(data, tracks(idx), f, ch,...
        'Handle', ha, 'iRange', dRange,...
        'Mode', ip.Results.Mode, 'DisplayType', ip.Results.DisplayType,...
        'ShowEvents', ip.Results.ShowEvents,...
        'ShowGaps', ip.Results.ShowGaps, 'Detection', detection, 'Colormap', cmap);
    
    axis(ha, 'off');
    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], [fpath 'frame' num2str(f, fmt) ext]);
    %print(h, '-djpeg100', '-loose', ['-r' num2str(zoom*72)], [fpath 'frame' num2str(f, fmt) ext]);
    cla(ha);
    fprintf('\b\b\b\b%3d%%', round(100*f/(nf)));
end
fprintf('\n');
close(h);

% Generate movie
fprintf('Generating movie... ');
fr = num2str(ip.Results.FrameRate);
cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'frame' fmt ext ' -b 20000k ' mpath ip.Results.FileName '.mp4'];
system(cmd);
fprintf('done.\n');
