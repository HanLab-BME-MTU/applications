function testScript_assaf(filename,params)

p = inputParser;
p.addRequired('filename', @(x)validateattributes(x,{'char'},{'nonempty'}));
p.addOptional('params', @(x)validateattributes(x,{'struct'},{'nonempty'}));

if nargin < 2
    params.timePerFrame = 5;
    params.nRois = 1;
    params.always = false;
    params.pixelSize = 1;
    params.patchSize = 15.0; % um
    params.nTime = floor(200 / params.timePerFrame); % frames to process (200 minutes)
    params.maxSpeed = 90; % um / hr (max cell speed)
    
    % for kymographs display
    params.kymoResolution.maxDistMu = 180; % um
    params.kymoResolution.min = params.patchSize;
    params.kymoResolution.stripSize = params.patchSize;
    params.kymoResolution.max = ceil(params.kymoResolution.maxDistMu/params.pixelSize);     
end

processTimeLapse(filename, params);
close all;
end