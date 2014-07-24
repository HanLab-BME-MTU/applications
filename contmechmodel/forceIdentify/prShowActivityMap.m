function prShowActivityMap(inData,varargin)
%prShowActivityMap: Display color coded activity map.
%
% SYNOPSIS: 
%    prShowActivityMap(inData);
%    prShowActivityMap(inData,'par1',value1);
%
% INPUT:
%    inData : Input data. It is expected to a 2D matrix (or a file that stores this matix) 
%             whose 1st dimension is space and 2nd dimension is time.
%
%    OPTIONAL PAR/VALUE pairs:
%      PAR                VALUE
%    --------       ------------------
%    'outDir'     : The path for storing the activity map. Default is the current path.
%    'outFileName': The name of output file. Default is the name of the variables stored in
%                   'inData' if it is a fiel. Otherwise, it is empty. When it is empty, the map will
%                   not be saved.
%    'cMap'       : Color map. Default, 'jet'.
%    'maxScore'   : The maximum score in the color range. Default, the maximum of data.
%    'minScore'   : The minimum score in the color range. Default, the minimum of data.
%    'divScore'   : A score that divide the score range into two parts
%                   whose color map are easily distinguished. This is
%                   useful when one tries to identify poly/depoly or
%                   protrusion/retraction events. Default, NaN (no
%                   division).
%
% OUTPUT:
%    The figure of the activity map is saved as ['map_' outFileName] in 'outDir'.
%
% AUTHOR:Lin Ji.
% DATE  : Sept. 18, 2006

figH = figure;

%Default parameters:
outDir   = pwd;
cMap     = colormap(jet(128));
maxScore = Inf;
minScore = -Inf;
divScore = NaN;

outFileName = '';

if nargin > 1
   for kk = 1:2:nargin-1
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
         case 'divScore'
            divScore = varargin{kk+1};
      end
   end
end

%Load the data.
if ischar(inData)
   s = load(inData);
   sFieldNames = fieldnames(s);
   A = s.(sFieldNames{1});
else
   A = inData;
end

if isinf(maxScore)
   maxScore = max(A(:));
end
if isinf(minScore)
   minScore = min(A(:));
end
%Repeat the last row and column since they will not be shown in 'pcolor'.
A = [A A(:,end)];
A = [A; A(end,:)];
A(1,end) = minScore;
A(2,end) = maxScore;
A(find(A>maxScore)) = maxScore;
A(find(A<minScore)) = minScore;

if ~isnan(divScore)
   if divScore >= maxScore || divScore <= minScore
      divScore = (maxScore+minScore)/2;
   end
   
   %Design the color map to separate the two score ranges. We use the lower
   %1/2 for small score range and upper 1/2 for big score range.
   numColors = size(cMap,1);
   startBigColor = ceil(numColors/2);
   endSmallColor = startBigColor-1;
   bigCMap = cMap(startBigColor:end,:);
   numBigColors   = size(bigCMap,1);
   numSmallColors = numBigColors*(divScore-minScore)/(maxScore-divScore);
   
   cSP = spapi(2,1:endSmallColor,cMap(1:endSmallColor,:).');
   smallCMap = fnval(cSP,linspace(1,endSmallColor,numSmallColors));
   divCMap = [smallCMap.'; bigCMap];
%    divCMap = [[linspace(cMap(1,1),cMap(endSmallColor,1),numSmallColors).' ...
%       linspace(cMap(1,2),cMap(endSmallColor,2),numSmallColors).' ...
%       linspace(cMap(1,3),cMap(endSmallColor,3),numSmallColors).']; divCMap];
else
   divCMap = cMap;
end

%pcolor(A(:,1:end-1));
figure(figH); hold off;
pcolor(A);
xlabel('Time step');
ylabel('Boundary segment');
title(['Activity map of ' outFileName]);
%axis([1 size(A,2) 1 size(A,1)+1]);
colormap(divCMap); colorbar;

if ~isempty(outFileName)
   mapFileName = [outDir filesep 'map_' outFileName '.fig'];
   saveas(figH,mapFileName,'fig');
end

