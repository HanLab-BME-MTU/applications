function [dispArrayFiles,imgIndexOfDTimePts] = getDispArrayFiles(dArrayDir)
% getDispArrayFiles: Get all tracked displacement array data files ('disp_array_###.mat')
%                    in 'dArrayDir'.
%
% INPUT:
%    dArrayDir: The directory that stores 'disp_array_###.mat'.
%
% OUTPUT:
%    dispArrayFiles    : Cell array of 'disp_array_###.mat' data files.
%    imgIndexOfDTimePts: Image index of dynamic time points given by the 
%                        index in the file name (sorted).
%
% AUTHOR: Lin Ji.
% DATE  : June 7, 2006

if nargout > 3
   error('Too many output arguments.');
end

%First, search for 'flowTrack' data file in the 'corr' directory of this
% project.
fileList  = dir([dArrayDir filesep '*.mat']);
fileNames = {fileList.name};
dArrayFile = fileNames(strmatch('disp_array',fileNames));

if length(dArrayFile) == 1
   imgIndexOfDTimePts = 1;
   dispArrayFiles     = dArrayFile;
   return;
end

imgIndexOfDTimePts = [];
%Find disp_array files whose names have the correct format:
% 'disp_array_###.mat'.
indPat = '_(\d+).mat';
vInd = [];
for kk = 1:length(dArrayFile)
   indToken = regexp(dArrayFile,indPat,'tokens');
   if length(indToken) == 1
      if ~isempty(indToken{1}{1})
         vInd = [vInd kk];
         imgIndexOfDTimePts = [imgIndexOfDTimePts str2num(indToken{1}{1})];
      end
   end
end

if ~isempty(vInd)
   [imgIndexOfDTimePts, sortI] = sort(imgIndexOfDTimePts);
   dispArrayFiles = dArrayFile(vInd(sortI));
else
   dispArrayFiles = {};
end

