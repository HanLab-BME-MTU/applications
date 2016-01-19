function [Pr,tform]=interplaneDeformation(MD,varargin)
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
p=ip.Results;

outputDir=[MD.outputDirectory_ '/stiching/'];mkdir(outputDir);

Pr=ExternalProcess(MD);
Pr.setInFilePaths({});
Pr.setOutFilePaths({[outputDir filesep 'stichingParameter.mat']});
pa = Pr.getParameters();
pa.name = 'interplaneDeformation';
pa.parameters = p;
Pr.setParameters(pa);
Pr.setDateTime();
MD.addProcess(Pr);

%% register channel, using first time point

if(ip.Results.computeShift)
    volFixed=MD.getChannel(1).loadStack(1);
    volReg=MD.getChannel(2).loadStack(1);
%%    
    [optim,metric]=imregconfig('monomodal');
    
    optim.MaximumStepLength=10;
    optim.MinimumStepLength=0.0001;
    [tform]=imregtform(volReg,volFixed,'translation',optim,metric);
    save([outputDir filesep 'stichingParameter.mat'],'tform');
end

%  %% Register 
% 
%  if(ip.Results.warp)
%      load([outputDir filesep 'driftParameter.mat']);
%      fileDir=[outputDir filesep 'registeredVol'];mkdir(fileDir);
%      parfor i=1:MD.nFrames_
%          jIdx=find((jumpIdx<i));
%          vol=MD.getChannel(1).loadStack(i);
%          for j=jIdx
%              vol=imwarp(vol,displacements{j},'nearest','OutputView',imref3d(size(vol)));
%          end
%          stackWrite(vol,[fileDir filesep 'registered-frame-'  num2str(i,'%04d') '.tif'])
%      end
%  end
