function nucleiStruc = singleNucleusSpotDetection(nucleiStruc, dataProperties, imageData, varargin)
% singleNucleusSpotDetection Detects spots from each multi-channel 3D stack
% of single nucleus
%   Detailed explanation goes here

% 05/2016 Ning Zhang

p = inputParser;
p.addRequired('nucleiStruc', @(x) ~isempty(x));
p.addRequired('channels', @(x) numel(x) > 1);
p.addRequired('imageData', @(x) ~isempty(x));
p.addParameter('detectionMethod', 'mnp', @isstr);
p.addParameter('mannualAdjMode', 0, @isnumeric);

p.parse(nucleiStruc, dataProperties.channel, imageData, varargin{:});

nucleiStruc = p.Results.nucleiStruc;
chaParams = p.Results.channels;
imageData = p.Results.imageData;
detectionMethod = p.Results.detectionMethod;
mannualAdjMode = p.Results.mannualAdjMode;

for chaNum = 1:numel(chaParams)
    chaName = chaParams(chaNum).name;
    switch lower(chaName)
        case 'dapi'
            continue;
            
        case 'green' 
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).green;
                spots = spotFind3D(nucStack, chaParams(chaNum), detectionMethod, mannualAdjMode);
                nucleiStruc(nucNum).greenSpot = spots.sp;
                clear nucStack
            end

        case 'red'
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).red;
                spots = spotFind3D(nucStack, chaParams(chaNum), detectionMethod, mannualAdjMode);
                nucleiStruc(nucNum).redSpot = spots.sp;
                clear nucStack
            end
            
        otherwise
            error('Unknown channels detected')
    end
   
end

spotsPlot3(nucleiStruc, imageData, dataProperties);