function compareDetection(MD,frameIdxs,xps,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('frameIdxs',@isnumeric);
ip.addRequired('xps')
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('printLabel',false,@islogical);
ip.addParamValue('showFrame',0,@isnumeric);
ip.parse(MD,frameIdxs,xps, varargin{:});

detectionSets=cell(1,numel(xps));
labels=cell(1,numel(xps));
for i=1:numel(xps)
    strParam=xps{i}{2}; strParam(cellfun(@isnumeric, strParam)) =  cellfun(@(x) sprintf('_%.4f', x), strParam(cellfun(@isnumeric, strParam)), 'UniformOutput', false);
    strParam=cell2mat(strParam);
    outputDir=[MD.outputDirectory_ filesep xps{i}{1} filesep strParam ];
    [MI,lab]=detectEB3(MD,'type',xps{i}{1},xps{i}{2}{:},'processFrames',frameIdxs,varargin{:});
    labels{i}=lab;
    mkdir(outputDir);
    if ip.Results.printLabel
        mkdir([outputDir filesep 'mask']);
        parfor tidx=frameIdxs
            stackWrite(lab{tidx},[outputDir filesep 'mask' filesep 'detect_T_' num2str(tidx,'%05d') '.tif']);
        end 
    end 
    save([outputDir filesep 'movieInfo.mat'],'MI');
    detectionSets{i}=MI;
end 

if(ip.Results.showFrame)
    vol=MD.getChannel(ip.Results.channel).loadStack(ip.Results.showFrame);
    imseriesmaskshow(vol,cellfun(@(x) x{ip.Results.showFrame},labels,'UniformOutput',false));
end