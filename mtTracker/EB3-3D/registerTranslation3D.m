function [Pr]=registerTranslation3D(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParamValue('channel',1,@isnumeric);
ip.addParamValue('computeImageDistance',true, @islogical);
ip.addParamValue('computeShift',true, @islogical);
ip.addParamValue('warp', true, @islogical);
ip.addParamValue('show', false, @islogical);
ip.parse(MD, varargin{:});

outputDir=[MD.outputDirectory_ '/regFile/'];mkdir(outputDir);

Pr=ExternalProcess(MD)
Pr.setInFilePaths({})
Pr.setOutFilePaths({[outputDir filesep 'driftParameter.mat']})
pa = Pr.getParameters();
pa.name = 'registration';
Pr.setParameters(pa);
Pr.setDateTime();
MD.addProcess(Pr);

%% detect a single jump
channelIdx=ip.Results.channel;
frameNb=MD.nFrames_;
dist=NaN(1,frameNb);
jumpIdx=0;
if(ip.Results.computeImageDistance)
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
     fileDir=[outputDir filesep 'registeredVol'];mkdir(fileDir);
     parfor i=1:MD.nFrames_
         jIdx=find((jumpIdx<i));
         vol=MD.getChannel(1).loadStack(i);
         for j=jIdx
             vol=imwarp(vol,displacements{j},'nearest','OutputView',imref3d(size(vol)));
         end
         stackWrite(vol,[fileDir filesep 'registered-frame-'  num2str(i,'%04d') '.tif'])
     end
 end