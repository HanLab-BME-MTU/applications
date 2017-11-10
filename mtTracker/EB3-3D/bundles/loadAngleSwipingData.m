function [fiberCell,randFiberCell]=loadAngleSwipingData(MD,varargin)
% Swipe different angle used to consider capute Detection, output the
% bundle, save them if a process is input. 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('process',[]);
ip.parse(MD,varargin{:});
p=ip.Results;

angleRange=0.02:0.01:0.1;

fiberCell=cell(1,length(angleRange));
randFiberCell=cell(1,length(angleRange));

parfor angleDistIdx=1:length(angleRange)

angleDist=angleRange(angleDistIdx);

%% Estimate bundle in KinPole axis
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles' filesep ['inlier_angle_' num2str(angleDist)]];
tmp=load([outputDirBundle filesep 'kin-MT-bundle.mat']);
fiberCell{angleDistIdx}=tmp.kinTracks;

outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles' filesep ['randInlier_angle_' num2str(angleDist)]];
tmp=load([outputDirBundle filesep 'kin-MT-bundle.mat']);
randFiberCell{angleDistIdx}=tmp.kinTracks;
end

displayAngleSwipeResults(fiberCell,randFiberCell)
