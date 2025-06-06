% makeMovieCME(data, varargin)
%
% Inputs: 
%             data : data structure from 'loadConditionData()'
%      {trackinfo} : track structure returned by loadtracks()
%
% Options: 
%       'Scalebar' : length of the scalebar. Default: 5 �m
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
        if iscell(data.framePaths{c})
            firstFrame = double(imread(data.framePaths{c}{1}));
            lastFrame = double(imread(data.framePaths{c}{data.movieLength}));
        else
            firstFrame = double(readtiff(data.framePaths{c}, 1));
            lastFrame = double(readtiff(data.framePaths{c}, data.movieLength));
        end
        dRange{c} = [min(min(firstFrame(:)),min(lastFrame(:))) max(max(firstFrame(:)),max(lastFrame(:)))];
    end
end

iptsetpref('ImshowBorder','tight');

% Figure options (PaperPosition is most important); points are not rendered !!
fopts = {'Visible', 'off', 'Position', [50 50 nx ny],...
    'InvertHardcopy', 'off', 'PaperUnits', 'Points', 'PaperSize', [nx ny],...
    'PaperPosition', [0 0 nx ny], 'PaperPositionMode', 'auto',...
    'DefaultLineLineSmoothing','on', 'DefaultPatchLineSmoothing','on'};

h = figure(fopts{:});

% Configure axes
ha = axes('Position', [0 0 1 1]);

frameopts = {'iRange', dRange,...
        'Mode', ip.Results.Mode, 'DisplayType', ip.Results.DisplayType,...
        'ShowEvents', ip.Results.ShowEvents,...
        'ShowGaps', ip.Results.ShowGaps};

% Generate frames
fmt = ['%0' num2str(ceil(log10(data.movieLength))) 'd'];
nf = numel(ip.Results.FrameRange);
fprintf('Generating movie frames:     ');
for f = ip.Results.FrameRange
    
    iopts = frameopts;
    
    if ~isempty(ip.Results.Detection)
        iopts = [iopts {'Detection', ip.Results.Detection(f)}];
    end
    
    if ~isempty(tracks)
        idx = [tracks.start]<=f & f<=[tracks.end];
        itracks = tracks(idx);
        if ~isempty(ip.Results.Colormap)
            iopts = [iopts {'ColorMap', ip.Results.Colormap(idx,:)}]; %#ok<*AGROW>
        end
    else
        itracks = [];
    end
    
    plotFrame(data, itracks, f, ch, 'Handle', ha, iopts{:});
    
    axis(ha, 'off');
    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], [fpath 'frame' num2str(f, fmt) ext]);
    %print(h, '-djpeg100', '-loose', ['-r' num2str(zoom*72)], [fpath 'frame' num2str(f, fmt) ext]);
    cla(ha);
    fprintf('\b\b\b\b%3d%%', round(100*f/(nf)));
end
fprintf('\n');
fprintf(['Frames saved to ' getShortPath(data) filesep 'Movies' filesep 'Frames.\n']);
close(h);

% Generate movie, if on a unix system with ffmpeg
if isunix && ~system('which ffmpeg >/dev/null 2>&1')
    fprintf('Generating movie ... ');
    fr = num2str(ip.Results.FrameRate);
    nx = 2*floor(nx/2);
    ny = 2*floor(ny/2);
    cmd = ['ffmpeg -y -r ' fr ' -i ' fpath 'frame' fmt ext ' -vf "scale=' num2str(nx) ':' num2str(ny)...
            '" -c:v libx264 -crf 22 -pix_fmt yuv420p ' mpath ip.Results.FileName '.mp4'];
    system(cmd);
    fprintf(' done.\n');
else
    fprintf('A unix system with ffmpeg installed is required to generate movies automatically.\n');
end
