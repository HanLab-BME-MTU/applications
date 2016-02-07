function [ output_args ] = GCAAddSplineFitFilopodiaMovie(movieData,varargin)
%GCAAddActinContentMetricMovie(movieData) 

%% CheckInput
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' ... 
    filesep 'StepsToReconstruct' filesep... 
    'VII_filopodiaBranch_fits' filesep 'Channel_1']; 
ip.addParameter('InputDirectory',defaultInDir,@(x) ischar(x)); 


ip.parse(varargin{:});
%%

load([ip.Results.InputDirectory filesep 'filoBranch.mat']); 

for iFrame = 1:length(filoBranch)-1
%for iFrame = 1:35
% extract the veil
%veilMask = analInfo(iFrame).masks.neuriteEdge;
% extract the img to feed into the function
%img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}])); 
% extract the filo info to read into the function 
filoInfo = filoBranch(iFrame).filoInfo;
% add the metric to the filo info - NOTE in the future might want to just
% calculate automaticaly at the time of fitting to be more efficient. 
filoInfo = GCASplineFitFilopodia(filoInfo); 
filoBranch(iFrame).filoInfo = filoInfo; 
%display(['Finished Frame ' num2str(iFrame)]); 

end
% resave the values 
save([ip.Results.InputDirectory ...
    filesep 'filoBranch.mat'],'filoBranch','-v7.3'); 





end

