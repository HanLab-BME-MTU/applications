function [Pr]=planeStitching(MD,registrationProcess,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('registrationProcess',@(x) isa(x,'ExternalProcess'));
%ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('computeImageDistance',true, @islogical);
ip.addParamValue('computeShift',true, @islogical);
ip.addParamValue('warp', true, @islogical);
ip.addParamValue('show', false, @islogical);
ip.parse(MD,registrationProcess, varargin{:});
p=ip.Results;

outputDir=[MD.outputDirectory_ '/stiching/'];mkdir(outputDir);

Pr=ExternalProcess(MD);
Pr.setInFilePaths({});
Pr.setOutFilePaths({[outputDir filesep 'stitchedMovieData.mat']});
pa = Pr.getParameters();
pa.name = 'planeStitching';
pa.parameters = p;
Pr.setParameters(pa);
Pr.setDateTime();
MD.addProcess(Pr);

 % Register 
 if(ip.Results.warp)
     registrationInfo=load(registrationProcess.outFilePaths_{1});
     fileDir=[outputDir filesep 'stitchedMovie'];mkdir(fileDir);
     mkdir([fileDir filesep  'ch1' filesep]);
     mkdir([fileDir filesep  'ch2' filesep]);
     mkdir([fileDir filesep  'analysis' filesep]);
     for i=1:MD.nFrames_
         vol1=MD.getChannel(1).loadStack(i);
         vol1=padarray(vol1,[1 size(vol1,2) size(vol1,3)]);
         vol2=MD.getChannel(2).loadStack(i);
         outputRef=imref3d(size(vol1));
%          outputRef.XWorldLimits(2) = outputRef.XWorldLimits(2)+500;
%          outputRef.YWorldLimits(2) = outputRef.YWorldLimits(2)+500;
%          outputRef.ZWorldLimits(2) = outputRef.ZWorldLimits(2)+1000;
         vol2=imwarp(vol2,registrationInfo.tform,'linear','OutputView',outputRef);
         stackWrite(vol1,[fileDir filesep  'ch1' filesep 'stichedChannel1'  num2str(i,'%04d') '.tif'])
         stackWrite(vol2,[fileDir filesep  'ch2' filesep 'stichedChannel2'  num2str(i,'%04d') '.tif'])
     end
     MD=MovieData();
     
     channelList=[Channel([fileDir  filesep 'ch1']),Channel([fileDir  filesep 'ch2'])];
     tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
     MD=MovieData(channelList,[fileDir filesep 'analysis'],'movieDataFileName_','stitchedMovieData.mat','movieDataPath_',[fileDir filesep 'analysis'], ...
         'pixelSize_',MD.pixelSize_,'pixelSizeZ_',MD.pixelSizeZ_,'timeInterval_',MD.timeInterval_);
     MD.setReader(tiffReader);
     MD.sanityCheck();
     MD.save();
 end
