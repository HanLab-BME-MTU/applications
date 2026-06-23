function filopodiaSegSweep(MD, varargin)
%FILOPODIASEGSWEEP  Diagnostic sweep for P1 body segmentation. On a single
%frame, re-segment the cell body under several threshold / morphology
%settings and overlay each resulting body outline on the talin image, in a
%grid, so you can pick the settings that capture the true cell edge.
%
%   filopodiaSegSweep(MD)
%   filopodiaSegSweep(MD, 'Frame',56, 'ChannelIndex',1)
%   filopodiaSegSweep(MD, 'Threshes',{'otsu','rosin'}, ...
%        'Scales',[1 0.6 0.5 0.4], 'Blurs',[1 2], ...
%        'Opens',[3 6 8], 'Closes',[4 8])
%
% By default it sweeps a threshold-method x scale grid at fixed morphology,
% which is the knob that matters most for the saturated-core talin images.
% Add 'Mode','morph' to instead sweep blur/open/close at a fixed threshold.
%
% Nothing is saved to the MovieData; this only writes a PNG/fig next to the
% current folder (or to 'OutDir' if given) for you to inspect.
% Sangyoon J. Han / 2026

ip = inputParser;
ip.addParameter('Frame',[]);
ip.addParameter('ChannelIndex',1,@isnumeric);
ip.addParameter('Mode','thresh',@ischar);    % 'thresh' | 'morph'
ip.addParameter('Threshes',{'otsu','rosin'},@iscell);
ip.addParameter('Scales',[1 0.6 0.5 0.4],@isnumeric);
ip.addParameter('Blurs',2,@isnumeric);
ip.addParameter('Opens',8,@isnumeric);
ip.addParameter('Closes',8,@isnumeric);
ip.addParameter('BodyMinArea',500,@isnumeric);
ip.addParameter('OutDir','',@ischar);
ip.parse(varargin{:});
R = ip.Results;

iChan = R.ChannelIndex;
fr = R.Frame; if isempty(fr), fr = max(1,round(MD.nFrames_/2)); end
img = double(MD.channels_(iChan).loadImage(fr));

% build the list of parameter sets to try
combos = struct('label',{},'thresh',{},'scale',{},'blur',{},'open',{},'close',{});
if strcmpi(R.Mode,'morph')
    th = R.Threshes{1}; sc = R.Scales(1);
    for b = R.Blurs
        for o = R.Opens
            for c = R.Closes
                combos(end+1) = mk(sprintf('%s blur%g op%g cl%g',th,b,o,c), ...
                    th, sc, b, o, c); %#ok<AGROW>
            end
        end
    end
else  % thresh mode
    b = R.Blurs(1); o = R.Opens(1); c = R.Closes(1);
    for ti = 1:numel(R.Threshes)
        for s = R.Scales
            lbl = sprintf('%s x%.2g', R.Threshes{ti}, s);
            combos(end+1) = mk(lbl, R.Threshes{ti}, s, b, o, c); %#ok<AGROW>
        end
    end
end

n = numel(combos);
nc = ceil(sqrt(n)); nr = ceil(n/nc);
h = figure('Color','w','Position',[80 80 360*nc 360*nr]);
for i = 1:n
    bm = segBody(img, combos(i), R.BodyMinArea);
    ax = subplot(nr,nc,i);
    imshow(img,[],'Parent',ax); hold(ax,'on'); colormap(ax,gray);
    B = bwboundaries(bm);
    for k = 1:numel(B)
        plot(ax,B{k}(:,2),B{k}(:,1),'-','Color',[1 0.2 0.2],'LineWidth',1.2);
    end
    hold(ax,'off');
    title(ax, combos(i).label, 'Interpreter','none','FontSize',9);
end
sgtitle(sprintf('P1 body segmentation sweep  (frame %d, ch %d)', fr, iChan));

outDir = R.OutDir; if isempty(outDir), outDir = pwd; end
if exist(outDir,'dir')~=7, mkdir(outDir); end
fbase = fullfile(outDir, sprintf('segSweep_%s_frame%04d', R.Mode, fr));
print(h,[fbase '.png'],'-dpng','-r150');
try, hgsave(h,fbase,'-v7.3'); catch, end
fprintf('Saved sweep -> %s.png\n', fbase);
end

% =====================================================================
function s = mk(label, thresh, scale, blur, op, cl)
s = struct('label',label,'thresh',thresh,'scale',scale, ...
    'blur',blur,'open',op,'close',cl);
end

% =====================================================================
function bodyMask = segBody(img, c, minArea)
% Replicates the body-segmentation block of segmentMovieFilopodia, with the
% threshold method/scale and morphology taken from combo c.
imgB = imgaussfilt(img, c.blur);
switch lower(c.thresh)
    case 'otsu',  level = thresholdOtsu(imgB);
    case 'rosin', level = thresholdRosin(imgB);
    otherwise,    level = str2double(c.thresh);
end
level = level * c.scale;
bodyMask = imgB > level;
bodyMask = imfill(bodyMask,'holes');
if any(bodyMask(:)), bodyMask = bwareafilt(bodyMask,1); end
if c.open  > 0, bodyMask = imopen(bodyMask, strel('disk',round(c.open))); end
if c.close > 0, bodyMask = imclose(bodyMask, strel('disk',round(c.close))); end
bodyMask = imfill(bodyMask,'holes');
bodyMask = bwareaopen(bodyMask, minArea);
end
