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
ip.addParamValue('mask',[],@isnumeric);
ip.addParamValue('jumpIdx',[],@isnumeric);
ip.addParamValue('computeImageDistance',true, @islogical);
ip.addParamValue('computeShift',true, @islogical);
ip.addParamValue('warp', true, @islogical);
ip.addParamValue('show', false, @islogical);
ip.addParamValue('warpMode','nearest', @ischar);
ip.parse(MD, varargin{:});
p=ip.Results;
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

% find mask offset (WARNING works only for cubic mask)
[maskMinX,maskMinY,maskMinZ]=ind2sub(size(p.mask), find(p.mask,1));
[maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(p.mask), find(p.mask,1,'last'));

%% detect a single jump
channelIdx=ip.Results.channel;
frameNb=MD.nFrames_;
dist=NaN(1,frameNb);
jumpIdx=ip.Results.jumpIdx;
if( (isempty(jumpIdx))||(ip.Results.computeImageDistance))
    parfor i=1:(frameNb-1)
        disp(['Processing frame ' int2str(i)]);
        voli=MD.getChannel(channelIdx).loadStack(i);
        volip1=MD.getChannel(channelIdx).loadStack(i+1);
        if(~isempty(p.mask))
            tmp=nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
            tmp(:)=voli(p.mask>0);
            voli=tmp;
            tmp(:)=volip1(p.mask>0);
            volip1=tmp;
        end          
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
    disp('Computing shift');
    displacements=cell(1,length(jumpIdx));
    fprintf('Frame registered :\n');
    parfor i=1:length(jumpIdx)       
        idx=jumpIdx(i)
        volFixed=MD.getChannel(channelIdx).loadStack(idx);
        volReg=MD.getChannel(channelIdx).loadStack(idx+1);
        if(~isempty(p.mask))
            tmp=nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
            tmp(:)=volFixed(p.mask>0);
            volFixed=tmp;
            tmp(:)=volReg(p.mask>0);
            volReg=tmp;
        end
        [optim,metric]=imregconfig('monomodal');
        [tform]=imregtform(volReg,volFixed,'translation',optim,metric);
        displacements{i}=tform;
        fprintf('.');
        if(mod(i,10)==0)
            fprintf('\n');
        end
    end
    save([outputDir filesep 'driftParameter.mat'],'jumpIdx','displacements');
end

%% Register 
 if(ip.Results.warp)
     load([outputDir filesep 'driftParameter.mat']);
     channelList=[];     
     for ch=1:length(MD.channels_)
         fileDir=[outputDir filesep 'registeredVol' filesep ip.Results.warpMode filesep 'ch' num2str(ch)];mkdir(fileDir);
         parfor i=1:MD.nFrames_
             jIdx=find((jumpIdx<i));
             vol=MD.getChannel(ch).loadStack(i);
             if(~isempty(p.mask))
                 tmp=nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
                 tmp(:)=vol(p.mask>0);
                 vol=tmp;
             end
             for j=jIdx
                 vol=imwarp(vol,displacements{j},ip.Results.warpMode,'OutputView',imref3d(size(vol)));
             end
             stackWrite(vol,[fileDir filesep 'registered-frame-'  num2str(i,'%04d') '.tif'])
             channelList=[channelList Channel([fileDir filesep])];              
         end
     end
     MD1=MovieData(channelList,[outputDir filesep 'analysis/'],'movieDataPath_' , outputDir, 'movieDataFileName_', 'movieData.mat');
     MD1.sanityCheck();
 end
 
 
function mkdir(path)
system(['mkdir -p ' path]);



