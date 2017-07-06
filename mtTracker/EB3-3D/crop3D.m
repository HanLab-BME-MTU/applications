function MDout=crop3D(MD,ROIorMask)    


if(numel(ROIorMask)>6)
    % find mask offset (WARNING works only for cubic mask)
    [maskMinX,maskMinY,maskMinZ]=ind2sub(size(ROIorMask), find(ROIorMask,1));
    [maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(ROIorMask), find(ROIorMask,1,'last'));
    ROIorMask=[maskMinX,maskMinY,maskMinZ,maskMaxX,maskMaxY,maskMaxZ];
end
ROI=ROIorMask;

outputDir=[MD.outputDirectory_ filesep 'cropped'];
MD.addProcess(ExternalProcess(MD,'crop'))
MD.getProcess(1).setInFilePaths({})
MD.getProcess(1).setOutFilePaths({outputDir});
p = MD.getProcess(1).getParameters();
p.metadata1 = 'crop3D';
p.metadata2 = ROI;
MD.getProcess(1).setParameters(p);
channelList=[];
for cIdx=1:length(MD.channels_)
    channelList=[channelList Channel([outputDir filesep 'ch' num2str(cIdx) filesep])];
    mkdirRobust([MD.getProcess(1).outFilePaths_{1} filesep 'ch' num2str(cIdx)]);
    parfor t=1:MD.nFrames_
        vol=MD.getChannel(cIdx).loadStack(t);
        stackWrite(vol(ROI(1):ROI(4),ROI(2):ROI(5),ROI(3):ROI(6)) ,[outputDir filesep 'ch' num2str(cIdx) filesep 'time-' num2str(t,'%04d') '.tif']);
    end 
end 

%% Make Movie Data
MDout=MovieData();
tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
MDout=MovieData(channelList,[outputDir filesep 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[outputDir filesep 'analysis'], ...
    'pixelSize_',MD.pixelSize_,'pixelSizeZ_',MD.pixelSizeZ_,'timeInterval_',MD.timeInterval_);
MDout.setReader(tiffReader);
MDout.sanityCheck();
MDout.save();

