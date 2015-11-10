function [Pr,jumpIdx,displacements]=registerTranslation3D(MD,varargin)
% Example of use with specified jumpIdx (whole stack here), and
% registration of the detection point stored in a movieInfo structure:
% =====
% [process,jumpIdx,displacements]=registerTranslation3D(MD,'channel',1,'jumpIdx',1:(MD.nFrames_-1));
% parfor i=1:MD.nFrames_
%     jIdx=find((jumpIdx<i));
%     for j=jIdx
%         movieInfo(i).xCoord(:,1)=movieInfo(i).xCoord(:,1)+displacements{j}.T(4,1);
%         movieInfo(i).yCoord(:,1)=movieInfo(i).yCoord(:,1)+displacements{j}.T(4,2);
%         movieInfo(i).zCoord(:,1)=movieInfo(i).zCoord(:,1)+displacements{j}.T(4,3);
%     end
% end
% save([outputDirDetect filesep 'detectionReg.mat'],'movieInfo');
% ======
% P. Roudot 2015

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('jumpIdx',[],@isnumeric);
ip.addParamValue('computeImageDistance',true, @islogical);
ip.addParamValue('computeShift',true, @islogical);
ip.addParamValue('warp', true, @islogical);
ip.addParamValue('show', false, @islogical);
ip.addParamValue('warpMode','nearest', @ischar);
ip.parse(MD, varargin{:});

outputDir=[MD.outputDirectory_ filesep 'regFile'];mkdir(outputDir);

Pr=ExternalProcess(MD)
Pr.setInFilePaths({})
Pr.setOutFilePaths({[outputDir filesep 'driftParameter.mat'],[outputDir filesep 'movieData.mat']})
pa = Pr.getParameters();
pa.name = 'registration';
pa.inputParam=ip.Results;
Pr.setParameters(pa);
Pr.setDateTime();
%MD.addProcess(Pr);

%% detect a single jump
channelIdx=ip.Results.channel;
frameNb=MD.nFrames_;
dist=NaN(1,frameNb);
jumpIdx=ip.Results.jumpIdx;
if( (~isempty(jumpIdx))||(ip.Results.computeImageDistance))
    parfor i=1:(frameNb-1)
        disp(['Processing frame ' int2str(i)]);
        voli=MD.getChannel(channelIdx).loadStack(i);
        volip1=MD.getChannel(channelIdx).loadStack(i+1);
        dist(i)=norm(double(volip1(:) - voli(:)));
    end
    thresh=nanmedian(dist)+3*1.4826*mad(dist,1);
    jumpIdx=find(dist>thresh)
    
    if(ip.Results.show)
        plot(dist,'color','r');
        hline(thresh);
    end
end

%% register image jumpIdx with jumpIdx+1
if(ip.Results.computeShift)
    displacements=cell(1,length(jumpIdx));
    parfor i=1:length(jumpIdx)
        idx=jumpIdx(i)
        volFixed=MD.getChannel(channelIdx).loadStack(idx);
        volReg=MD.getChannel(channelIdx).loadStack(idx+1);
        [optim,metric]=imregconfig('monomodal');
        [tform]=imregtform(volReg,volFixed,'translation',optim,metric);
        displacements{i}=tform;
    end
    save([outputDir filesep 'driftParameter.mat'],'jumpIdx','displacements');
end

%% Register 
 if(ip.Results.warp)
     load([outputDir filesep 'driftParameter.mat']);
     fileDir=[outputDir filesep 'registeredVol' filesep ip.Results.warpMode];mkdir(fileDir);
     parfor i=1:MD.nFrames_
         jIdx=find((jumpIdx<i));
         vol=MD.getChannel(1).loadStack(i);
         for j=jIdx
             vol=imwarp(vol,displacements{j},ip.Results.warpMode,'OutputView',imref3d(size(vol)));
         end
         stackWrite(vol,[fileDir filesep 'registered-frame-'  num2str(i,'%04d') '.tif'])
     end
     MD1=MovieData(Channel([fileDir filesep]),[outputDir filesep 'xp/'],'movieDataPath_' , outputDir, 'movieDataFileName_', 'movieData.mat')
     MD1.sanityCheck();
 end
