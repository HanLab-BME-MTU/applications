function [ output_args ] = GCAAddFilopodiaActinContentMetricMovie(movieData)
%GCAAddActinContentMetricMovie(movieData) 

load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' ...
    filesep 'analInfoTestSave.mat']); 




for iFrame = 1:length(analInfo)-1
% extract the veil
veilMask = analInfo(iFrame).masks.neuriteEdge;
% extract the img to feed into the function
img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}])); 
% extract the filo info to read into the function 
filoInfo = analInfo(iFrame).filoInfo;
% add the metric to the filo info - NOTE in the future might want to just
% calculate automaticaly at the time of fitting to be more efficient. 
[filoInfo,normFactPerFrame] = GCAAddFilopodiaActinContentMetric(img,veilMask,filoInfo); 
analInfo(iFrame).filoInfo = filoInfo; 
paramC{iFrame} = normFactPerFrame; 
end
% resave the values 
save([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1'...
    filesep 'analInfoTestSave.mat'],'analInfo','-v7.3'); 

%% save the normalization factor in the measurements folder 

% NOTE For final change to MEASUREMENT_EXTRACTION 
expFolder = ['PARAMETER_EXTRACTION' filespe 'Descriptor' filespe 'GrowthCone' filesep 'ExpressionNormalization']; 

if ~isdir(measurementFolder) 
    mkdir(measurementFolder);
end 
save([expFolder filesep 'param_ExpressionNormalization.mat'],'paramC'); 






end

