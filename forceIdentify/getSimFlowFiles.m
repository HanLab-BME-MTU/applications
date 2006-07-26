function [simFlowFiles,imgIndexOfDTimePts] = getSimFlowFiles(simFlowDir)
% getSimFlowFiles: Get all simulated flow data files ('iDispField###.mat')
%                    in 'simFlowDir'.
%
% INPUT:
%    simFlowDir: The directory that stores 'iDispField###.mat'.
%
% OUTPUT:
%    simFlowFiles    : Cell array of 'iDispField###.mat' data files.
%    imgIndexOfDTimePts: Image index of dynamic time points given by the 
%                        index in the file name (sorted).
%
% AUTHOR: Lin Ji.
% DATE  : July 19, 2006

if nargout > 3
   error('Too many output arguments.');
end

%First, search for 'flowTrack' data file in the 'corr' directory of this
% project.
fileList  = dir([simFlowDir filesep '*.mat']);
fileNames = {fileList.name};
sFlowFile = fileNames(strmatch('iDispField',fileNames));

if length(sFlowFile) == 1
   imgIndexOfDTimePts = 1;
   simFlowFiles     = sFlowFile;
   return;
end

imgIndexOfDTimePts = [];
%Find disp_array files whose names have the correct format:
% 'disp_array_###.mat'.
indPat = 'iDispField(\d+).mat';
vInd = [];
for kk = 1:length(sFlowFile)
   indToken = regexp(sFlowFile,indPat,'tokens');
   if length(indToken) == 1
      if ~isempty(indToken{1}{1})
         vInd = [vInd kk];
         imgIndexOfDTimePts = [imgIndexOfDTimePts str2num(indToken{1}{1})];
      end
   end
end

if ~isempty(vInd)
   [imgIndexOfDTimePts, sortI] = sort(imgIndexOfDTimePts);
   simFlowFiles = sFlowFile(vInd(sortI));
else
   simFlowFiles = {};
end

