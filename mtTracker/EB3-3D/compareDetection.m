function compareDetection(MD,frameIdx,xps,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('frameIdx',@isnumeric);
ip.addRequired('xps')
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('printLabel',false,@logical);
ip.parse(MD,frameIdx,xps, varargin{:});

detectionSetss=cell(1,numel(xps));
labels=cell(1,numel(xps));
for i=1:numel(xps)
    [detectionSets{i},lab]=detectEB3(MD,'type',xps{i}{1},xps{i}{2}{:},'processFrames',[frameIdx],varargin{:});
    labels{i}=lab>0;
    if ip.Results.printLabel
        stackWrite(labels{i},['detect' xps{i}{1} '.tif']);
    end 
end 

vol=MD.getChannel(ip.Results.channel).loadStack(frameIdx);
imseriesmaskshow(vol,labels)
