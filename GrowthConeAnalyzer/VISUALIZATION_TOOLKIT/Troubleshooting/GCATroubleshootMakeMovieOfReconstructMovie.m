function [ output_args ] = GCATroubleshootMakeMovieOfReconstructMovie(movieData,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% INPUTPARSER
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultInDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits'];

defaultOutDir = [movieData.outputDirectory_ ]; 

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('InputDirectory', defaultInDir,@(x) ischar(x)); 

ip.addParameter('frames',1)
ip.addParameter('outDirType','perFrame'); 


ip.parse(varargin{:});
%% Initiate 

imDir  = movieData.getChannelPaths{1}; 
frames = ip.Results.frames; 
%% Wrap 
for iFrame = 1:length(frames) 
    if strcmpi(ip.Results.outDirType,'perFrame');  
saveDir = [ip.Results.OutputDirectory filesep 'Reconstruct_Movies' filesep 'Frame_' num2str(frames(iFrame),'%03d')]; 
    else 
        saveDir = [ip.Results.OutputDirectory filesep 'Reconstruct_Movie']; 
    end 

if ~isdir(saveDir)
    mkdir(saveDir)
end 

load([movieData.outputDirectory_ filesep ...
    'SegmentationPackage' filesep 'StepsToReconstruct' ...
    filesep 'IV_veilStem_length' filesep 'Channel_1' filesep 'veilStem.mat']); 
%  load([movieData.outputDirectory_ filesep  'SegmentationPackage' filesep ... 
%     'StepsToReconstruct' filesep 'VI_filopodiaBranch_reconstruction' filesep 'Channel_1' filesep 'filoBranch.mat']); 

% load([movieData.outputDirectory_ filesep 'SegmentationPackage' filesep ...
%     'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep 'Channel_1' filesep 'filoBranch.mat']);

load([ip.Results.InputDirectory filesep 'filoBranch.mat'])


pixSize_um = movieData.pixelSize_/1000; 

%load([movieData.outputDirectory_ filesep 'filopodia_reconstruct' filesep 'Filopodia_Reconstruct_Channel_1' filesep 'analInfoTestSave.mat']); 
 
%load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' filesep 'analInfoTestSave.mat'])
[hSet filoFilterSet,filterParams] = GCATroubleShootMakeMovieOfReconstruct(filoBranch,veilStem,frames(iFrame),pixSize_um,imDir); 
filoInfo = filoBranch(iFrame).filoInfo; 
filoFilterSet = filoFilterSet{iFrame}; 

count = 1; 
for i = 1:length(hSet) 
    % double it up to slow it down 
    saveas(hSet(i).h, [saveDir filesep num2str(count,'%03d') '.png']);
    count= count+1; 
    saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']); 
    count = count+1; 
    
    % save two more of the last frame because windows media player clips
    % the last two frames 
    
    if i == length(hSet); 
       
        saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']); 
        saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']);
    end 
         
    
end 

save([saveDir filesep 'filoInfo.mat' ],'filoFilterSet','filterParams','filoInfo'); 


execute = ['ffmpeg -r 1 -i ' saveDir filesep '%03d.png' ...
    ' -b 2000k ' saveDir filesep 'ReconstructMovie' num2str(frames(iFrame),'%03d') '.wmv'];
system(execute);

execute = ['ffmpeg -r 1 -i ' saveDir filesep '%03d.png' ...
    ' -b 2000k ' saveDir filesep 'ReconstructMovie' num2str(frames(iFrame),'%03d') '.mp4'];
system(execute);
end 



end

