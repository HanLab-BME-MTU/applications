function [ output_args ] = GCAAddFilopodiaCurvatureMetricMovie(movieData)
%GCAAddActinContentMetricMovie(movieData) 

load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' ...
    filesep 'analInfoTestSave.mat']); 

for iFrame = 1:length(analInfo)-1
% extract the veil
%veilMask = analInfo(iFrame).masks.neuriteEdge;
% extract the img to feed into the function
%img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}])); 
% extract the filo info to read into the function 
filoInfo = analInfo(iFrame).filoInfo;
% add the metric to the filo info - NOTE in the future might want to just
% calculate automaticaly at the time of fitting to be more efficient. 
filoInfo = GCAAddFilopodiaCurvature(filoInfo); 
analInfo(iFrame).filoInfo = filoInfo; 

end
% resave the values 
save([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1'...
    filesep 'analInfoTestSave.mat'],'analInfo','-v7.3'); 





end

