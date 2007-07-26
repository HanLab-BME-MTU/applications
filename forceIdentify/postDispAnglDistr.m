function postDispAnglDistr(varargin)
%This function displays the distribution of angles between flow and force vectors in a circle.
%
% SYNOPSIS:
%    To use default number of bins for histogram:
%    postDispAnglDistr;
%
%    or,
%    postDispAnglDistr(numBins);
%

numBins = 40;
if nargin == 1
   numBins = varargin{1};
end

baseVar = evalin('base','who');
for kk = 1:length(baseVar)
   BV.(baseVar{kk}) = evalin('base',baseVar{kk});
end

if ~isdir([BV.reslDir filesep 'anglDistrFig'])
   success = mkdir(reslDir,'anglDistrFig');
   if ~success
      error('trouble making directory');
   end
end
anglDistrFigDir = [BV.reslDir filesep 'anglDistrFig'];

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:BV.numDTimePts;
else
   selTimeSteps = answer;
end

figH = figure;
backStr = '';
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   for kk = 1:length(backStr)
      fprintf(1,'\b');
   end
   backStr = sprintf('    Time step %d ...', jj);
   fprintf(1,backStr);

   imgIndex = BV.imgIndexOfDTimePts(jj);
   indexStr = sprintf(BV.imgIndexForm,imgIndex);

   forceFieldFile = [BV.forceFieldDir filesep 'forceField' ...
      indexStr '.mat'];
   s = load(forceFieldFile);
   forceField = s.forceField;

   anglDistrFigFile = [anglDistrFigDir filesep 'anglDistrFig' indexStr '.fig'];

   %Show histogram of angle distribution in 50 bins.
   titleStr = sprintf(['Distribution of angles between flow and force vectors\n' ...
      'Time Step: %d Image Index: %d '],jj,imgIndex);
   figure(figH);
   [binN binX] = hist(forceField.angleFlowBF*180/pi,numBins);
   bar(binX,binN/length(forceField.angleFlowBF),0.95);
   title(titleStr);
   saveas(figH,anglDistrFigFile,'fig');
end
fprintf(1,'\n');

