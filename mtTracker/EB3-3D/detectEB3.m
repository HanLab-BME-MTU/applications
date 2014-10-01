function [movieInfo,label] = detectEB3(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('processFrames',[], @isnumeric);
ip.addParamValue('waterThresh', 120, @isnumeric);
ip.addParamValue('waterStep', 10, @isnumeric);
ip.addParamValue('showAll', false, @islogical);
ip.addParamValue('type', 'watershed',  @ischar);
ip.parse(MD, varargin{:});


processFrames=[];
if isempty(ip.Results.processFrames)
    processFrames=1:numel(MD.getChannel(ip.Results.channel).getImageFileNames);
else
    processFrames=ip.Results.processFrames
end

for timePoint=processFrames
    vol=MD.getChannel(ip.Results.channel).loadStack(timePoint);

    switch ip.Results.type
        case 'watershed'
           [movieInfo,label]=detectComets3D(vol,10,ip.Results.waterThresh,[1 1 1]);
        otherwise 
           disp('Unsupported detection method.');
           disp('Supported method:');
           disp('\twatershed');
    end
    
    if ip.Results.showAll
        stackShow(vol,'overlay',label);
        imseriesmaskshow(vol,label);
    end
end 

% $$$ for frameIndex=ip.Results.processFrames 
% $$$ 
% $$$ if ~strcmp(ip.Results.type,'watershed')
% $$$     [movieInfo,label]=detectComets3D(vol,10,120,[1 1 1]);
% $$$ end 