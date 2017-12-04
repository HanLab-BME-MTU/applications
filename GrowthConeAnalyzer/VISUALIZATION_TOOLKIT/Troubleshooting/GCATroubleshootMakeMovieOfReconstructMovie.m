function [ output_args ] = GCATroubleshootMakeMovieOfReconstructMovie(movieData,varargin)
% GCATroubleshootMakeMovieOfReconstructMovie


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
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep 'Channel_1'];

defaultInDirVeilStem = [movieData.outputDirectory_ filesep ...
    'SegmentationPackage' filesep 'StepsToReconstruct' ...
    filesep 'IV_veilStem_length' filesep 'Channel_1']; 

defaultOutDir = [movieData.outputDirectory_ ]; 

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('InputDirectoryFiloBranch', defaultInDir,@(x) ischar(x)); 
ip.addParameter('InputDirectoryVeilStem',defaultInDirVeilStem,@(x) ischar(x)); 

ip.addParameter('frames',1)
ip.addParameter('outDirType','perFrame'); 

%ip.addParameter('figFiles',true); % need to add 
ip.addParameter('epsFiles',false); 

ip.addParameter('makeMovies',false); % flag to make movies : Requires ffmpeg

ip.addParameter('writeTitles',true); 
ip.addParameter('screen2png',false); 
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


load( [ip.Results.InputDirectoryVeilStem filesep 'veilStem.mat']); 


load([ip.Results.InputDirectoryFiloBranch filesep 'filoBranch.mat']);

pixSize_um = movieData.pixelSize_/1000; 

inputFrames = frames(iFrame); 

[hSet filoFilterSet,filterParams] = GCATroubleShootMakeMovieOfReconstruct(filoBranch,veilStem,inputFrames,pixSize_um,imDir,...
  'writeTitles',ip.Results.writeTitles); 
    
%filoInfo = filoBranch(inputFrames).filoInfo; 

count = 1; 
for i = 1:length(hSet) 
    % double it up to slow it down 
%     saveas(hSet(i).h, [saveDir filesep num2str(count,'%03d') '.png']);
%     count= count+1;
    cFileName = [saveDir filesep num2str(count,'%03d') '.png'];

    if ip.Results.screen2png
        helperScreen2png(cFileName,'figureHandle',hSet(i).h);
    else
        saveas(hSet(i).h,cFileName);
    end
      if ip.Results.epsFiles 
      saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.eps'],'psc2'); 
      end 
    count = count+1; 
    
    % save two more of the last frame because windows media player clips
    % the last two frames 
%     
%     if i == length(hSet)-2; 
%        
%         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']); 
%         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']);
%         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.fig']); 
%         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.eps'],'psc2'); 
%     end 
  
end 
%% make the separate directory 
% extraDir = [saveDir filesep 'NiceFigs'];
% if ~isdir(extraDir) 
%     mkdir(extraDir) ; 
% end 
% 
% names{1} = 'ColorByBrewer'; 
% names{2} = 'ColorByJet'; 
% m = length(hSet)-1:length(hSet); 
% for i =1:2 
%      
%     if ip.Results.screen2png
%         helperScreen2png([extraDir filesep names{i} '.png'],'figureHandle',hSet(m(i)).h);
%     else
%         saveas(hSet(m(i)).h,[extraDir filesep names{i} '.png']);
%     end
%     
%     
%      saveas(hSet(m(i)).h,[extraDir filesep names{i} '.eps'],'psc2');
%      saveas(hSet(m(i)).h,[extraDir filesep names{i} '.fig']);   
% end     
%% 
% imgSize = movieData.imSize_; 

% save only the current filoInfo filterSet and filterParams 
% save([saveDir filesep 'filoInfo.mat' ],'filoFilterSetC','filterParams','filoInfoG','filoInfo','imgSize'); 

load([ip.Results.InputDirectoryFiloBranch filesep 'params.mat']); 
copyfile([ip.Results.InputDirectoryFiloBranch  filesep 'params.mat'],[saveDir filesep 'paramsRec.mat']); 

if ip.Results.makeMovies
%     execute = ['ffmpeg -r 1 -i ' saveDir filesep '%03d.png' ...
%         ' -b 2000k ' saveDir filesep 'ReconstructMovie' num2str(frames(iFrame),'%03d') '.wmv'];
%     system(execute);
    
    execute = ['ffmpeg -r 1 -i ' saveDir filesep '%03d.png' ...
        ' -b 2000k ' '-pix_fmt' ' yuv420p ' saveDir filesep 'ReconstructMovie' num2str(frames(iFrame),'%03d') '.mp4'];
    system(execute);
end
end 

end

