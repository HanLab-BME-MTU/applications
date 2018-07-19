%imxTimestamp() uses ImageMagick to add a time stamp to movie frames
%
% For supported 'color' values, see http://www.imagemagick.org/script/color.php 

% Francois Aguet, 08/2013

function imxTimestamp(paths, t, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('paths');
ip.addRequired('t');
ip.addParamValue('FontSize', 20);
ip.addParamValue('Location', 'SouthWest');
ip.addParamValue('Offsets', [100 20]);
ip.addParamValue('Crop', []);
ip.addParamValue('Color', 'white'); % see 
ip.parse(paths, t, varargin{:});
crop = ip.Results.Crop;

if ischar(paths)
    paths = {paths};
end

info = imfinfo(paths{1});
nx = info.Width;
ny = info.Height;

dx = ip.Results.Offsets(1);
dy = ip.Results.Offsets(2);
wx = 11;

switch ip.Results.Color
    case 'w'
        
    case 'k'
        
end

for i = 1:numel(paths)
    str = num2str(t(i), '%.0f');
    cmd = ['export DYLD_LIBRARY_PATH=""; convert ' paths{i} ...
        ' -fill black -gravity SouthEast -font Helvetica -pointsize '  num2str(ip.Results.FontSize)];

    if ~isempty(crop)
                cmd = [cmd ' -crop ' num2str(round(nx*crop(3))) 'x' num2str(round(ny*crop(4))) ...
            '+' num2str(round(nx*(1-crop(1)-crop(3)))) '+' num2str(round(ny*crop(2))) ' +repage '];
        nx = round(nx*crop(3));
        ny = round(ny*crop(4));
    end

    cmd = [cmd ' -fill ' ip.Results.Color ' -annotate +' num2str(nx-dx) '+' num2str(dy) ' "s" ']; %#ok<*AGROW>
    for k = numel(str):-1:1
        cmd = [cmd ' -fill ' ip.Results.Color ' -annotate +' num2str(nx-dx+wx*(k)) '+' num2str(dy) ' "' str(numel(str)-k+1) '" '];
    end
    cmd = [cmd paths{i}];
    system(cmd);
end
