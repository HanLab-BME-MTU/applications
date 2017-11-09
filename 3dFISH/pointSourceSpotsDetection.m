function nucleiStruc = pointSourceSpotsDetection(nucleiStruc, dataProperties, varargin)
% pointSourceSpotsDetection Detects spots from each multi-channel 3D stack
% of single nucleus using Francois method
%   Detailed explanation goes here

% 07/2016 Ning Zhang

p = inputParser;
p.addRequired('nucleiStruc', @(x) ~isempty(x));
p.addRequired('channels', @(x) numel(x) > 1);
p.parse(nucleiStruc, dataProperties.channel, varargin{:});

nucleiStruc = p.Results.nucleiStruc;
chaParams = p.Results.channels;

for chaNum = 1:numel(chaParams)
    chaName = chaParams(chaNum).name;
    switch lower(chaName)
        case 'dapi'
            continue;
            
        case 'green' 
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).green;                
                [pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(nucStack, [chaParams(chaNum).psfSigmaMD, chaParams(chaNum).psfSigma(2)]);
                for spotNum = 1:numel(pstruct.x)
                    nucleiStruc(nucNum).greenSpot(spotNum).cord = [pstruct.x(spotNum), pstruct.y(spotNum), pstruct.z(spotNum)];
                end
                clear nucStack
            end

        case 'red'
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).red;
                [pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(nucStack, [chaParams(chaNum).psfSigmaMD, chaParams(chaNum).psfSigma(2)]);
                for spotNum = 1:numel(pstruct.x)
                    nucleiStruc(nucNum).redSpot(spotNum).cord = [pstruct.x(spotNum), pstruct.y(spotNum), pstruct.z(spotNum)];
                end
                clear nucStack
            end
            
        otherwise
            error('Unknown channels detected')
    end
   
end