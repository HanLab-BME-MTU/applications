function [siCorrM,siShiftJ,siShiftI] = prShowCorrScoreMap(corrM,shiftJ,shiftI,varargin)
%prShowCorrScoreMap: Display color coded cross-correlation score map.
%
% SYNOPSIS: 
%    [siCorrM,siShiftJ,siShiftI] = prShowCorrScoreMap(corrM,shiftJ,shiftI);
%    [siCorrM,siShiftJ,siShiftI] = prShowCorrScoreMap(corrM,shiftJ,shiftI,'par1',value1,...);
%
% INPUT:
%    corrM: A 2D matrix that stores the cross-correlation score between two activity map.
%           The first dimension is shift in space and 2nd dimension is shift in time.
%    shiftJ: Shift (integer) in J-direction (vertical or Y-axis, 1st dimension of 'corrM').
%    shiftI: Shift (integer) in I-direction (horizontal or X-axis, 2nd dimension of 'corrM').
%
%    OPTIONAL PAR/VALUE pairs:
%      PAR                VALUE
%    --------       ------------------
%    'outDir'     : The path for storing the correlation score map. Default is the current path.
%    'outFileName': The name of output file. Default is the name of the variables stored in
%                   'inData' if it is a file. Otherwise, it is empty. 
%                   Note: when it is empty, the map will not be saved.
%    'cMap'       : Color map. Default, 'jet'.
%    'maxScore'   : The maximum score in the color range. Default, the maximum of data.
%    'minScore'   : The minimum score in the color range. Default, the minimum of data.
%    'numSubIntSamples': Number of sub-integer samples. Default, 4.
%    'showGridLines': 'yes' or 'no' (default).
%    'titleStr'   : Specify a title.
%
% OUTPUT:
%    siCorrM : The correlation score matrix resampled to sub-integer shift by B-spline 
%              interpolation.
%    siShiftJ: Shift (at sub-integer) in J-direction (vertical or Y-axis, 1st dimension of 'corrM').
%    siShiftI: Shift (at sub-integer) in I-direction (horizontal or X-axis, 2nd dimension of 'corrM').
%
%    The figure of the correlation map is shown and saved as ['corr_' outFileName '.fig'] 
%    in 'outDir'. The correlation score matrix resampled to sub-integer shift will also be 
%    saved as ['siCorr_' outFileName '.mat'] in 'outDir'.
%
% AUTHOR:Lin Ji.
% DATE  : Sept. 19, 2006

%Make sure the size of 'corrM' matches that of 'shiftJ' and 'shiftI'.
if size(corrM,1) ~= length(shiftJ) || size(corrM,2) ~= length(shiftI)
   error('The size of ''corrM'' does not match that of ''shiftJ'' and ''shiftI.''');
end

%Default parameters:
outDir           = pwd;
cMap             = 'jet';
maxScore         = Inf;
minScore         = -Inf;
numSubIntSamples = 4;
showGridLines    = 'no';
titleStr         = 'Correlation map';

outFileName = '';

if nargin > 3 
   for kk = 1:2:nargin-3
      switch varargin{kk}
         case 'outDir'
            outDir = varargin{kk+1};
         case 'outFileName'
            outFileName = varargin{kk+1};
         case 'cMap'
            cMap = varargin{kk+1};
         case 'maxScore'
            maxScore = varargin{kk+1};
         case 'minScore'
            minScore = varargin{kk+1};
         case 'numSubIntSamples'
            numSubIntSamples = varargin{kk+1};
         case 'showGridLines'
            showGridLines = varargin{kk+1};
         case 'titleStr'
            titleStr = varargin{kk+1};
      end
   end
end

if numSubIntSamples == 0
   siShiftJ = shiftJ;
   siShiftI = shiftI;
   siCorrM  = corrM;
else
   %Spline interpolation of integer shifts.
   sp = spapi({4 4},{shiftJ,shiftI},corrM);

   %Add sub-integer samples to 'shiftJ' and 'shiftI'.
   siShiftJ = [];
   for kk = 1:length(shiftJ)-1
      siShiftJ = [siShiftJ linspace(shiftJ(kk),shiftJ(kk+1),numSubIntSamples+2)];
      siShiftJ(end) = [];
   end
   siShiftJ = [siShiftJ shiftJ(end)];

   siShiftI = [];
   for kk = 1:length(shiftI)-1
      siShiftI = [siShiftI linspace(shiftI(kk),shiftI(kk+1),numSubIntSamples+2)];
      siShiftI(end) = [];
   end
   siShiftI = [siShiftI shiftI(end)];

   %Evaluate at sub-integer shifts.
   siCorrM = fnval(sp,{siShiftJ,siShiftI});
end

%Save the sub-integer sampled correlation matrix.
if ~isempty(outFileName)
   siCorrFileName = [outDir filesep 'siCorr_' outFileName '.mat'];
   save(siCorrFileName,'siCorrM','siShiftJ','siShiftI');
end

%Show the (sub-integer sampled) correlation score map.
if isinf(maxScore)
   maxScore = max(siCorrM(:));
end
if isinf(minScore)
   minScore = min(siCorrM(:));
end

%To make the color range to be [minScore maxScore], we extend 'siCorrM' to cut off out of range
%values.
siShiftI = [siShiftI 2*siShiftI(end)-siShiftI(end-1)];
siShiftJ = [siShiftJ 2*siShiftJ(end)-siShiftJ(end-1)];
siCorrM = [siCorrM siCorrM(:,end)];
siCorrM = [siCorrM; siCorrM(end,:)];
siCorrM(1,end) = minScore;
siCorrM(2,end) = maxScore;
siCorrM(find(siCorrM>maxScore)) = maxScore;
siCorrM(find(siCorrM<minScore)) = minScore;

figH = figure;
surfH = pcolor(siShiftI,siShiftJ,siCorrM);
if strcmp(showGridLines,'no')
   set(surfH,'LineStyle','none');
end
set(gca,'XTick',shiftI);
set(gca,'YTick',shiftJ);
xlabel('Time step shift');
ylabel('Boundary segment shift');
title(titleStr);
colormap(cMap); colorbar;
%axis([1 siShiftI(end-1) 1 siShiftJ(end)]);

if ~isempty(outFileName)
   mapFileName = [outDir filesep 'siCorr_' outFileName '.fig'];
   saveas(figH,mapFileName,'fig');
end

