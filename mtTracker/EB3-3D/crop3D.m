function MDout=crop3D(MD,ROIorMask,varargin)    
ip = inputParser; ip.CaseSensitive = false;  ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('ROIorMask',@isnumeric);
ip.addParameter('name', '',@ischar);
ip.parse(MD,ROIorMask, varargin{:});
p=ip.Results;

if(numel(ROIorMask)>6)
    % find mask offset (WARNING works only for cubic/rectangular mask)
    [maskMinX,maskMinY,maskMinZ]=ind2sub(size(ROIorMask), find(ROIorMask,1));
    [maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(ROIorMask), find(ROIorMask,1,'last'));
    maskMaxSize=min([maskMaxX,maskMaxY,maskMaxZ],[MD.imSize_,MD.zSize_]);
    maskMaxX=maskMaxSize(1); maskMaxY=maskMaxSize(2); maskMaxZ=maskMaxSize(3);
    ROIorMask=[maskMinX,maskMinY,maskMinZ,maskMaxX,maskMaxY,maskMaxZ];
end
ROI=ROIorMask;

outputDir=[MD.outputDirectory_ filesep 'cropped' filesep p.name filesep];
mkdirRobust(outputDir);
process=ExternalProcess(MD,'crop3D');
process.setInFilePaths({});
process.setOutFilePaths({[outputDir filesep 'analysis' filesep 'movieData.mat']});
process.setProcessTag(['crop3D_' p.name])
p = process.getParameters();
p.metadata1 = 'crop3D';
p.metadata2 = ROI;
process.setParameters(p);
MD.addProcess(process);
channelList=[];
for cIdx=1:length(MD.channels_)
    channelList=[channelList Channel([outputDir filesep 'ch' num2str(cIdx) filesep])];
    mkdirRobust([outputDir filesep 'ch' num2str(cIdx)]);
    parfor t=1:MD.nFrames_
        vol=MD.getChannel(cIdx).loadStack(t);
        stackWrite(vol(ROI(1):ROI(4),ROI(2):ROI(5),ROI(3):ROI(6)) ,[outputDir filesep 'ch' num2str(cIdx) filesep 'time-' num2str(t,'%04d') '.tif']);
    end 
end 

%% Make Movie Data
mkdirRobust([outputDir filesep 'analysis']);
MDout=MovieData();
tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
MDout=MovieData(channelList,[outputDir filesep 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[outputDir filesep 'analysis'], ...
    'pixelSize_',MD.pixelSize_,'pixelSizeZ_',MD.pixelSizeZ_,'timeInterval_',MD.timeInterval_);
MDout.setReader(tiffReader);
MDout.sanityCheck();
MDout.save();

