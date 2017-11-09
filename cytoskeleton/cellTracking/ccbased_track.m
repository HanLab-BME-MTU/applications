function [d,ccMap,x1,x2] = ccbased_track(tImg,tPos,tDim,sImg,dmax,subpix,d0)
%CCTRACK correlation-based feature tracking
%
% all vectors are given in left-handed coordinate system,
% i.e. ----> x(1)
%      |
%      |
%      \/ x(2)
%
% SYNOPSIS [d,ccMap] = cctrack(tImg,tPos,tDim,sImg,dmax,subpix,d0)
% 
% INPUT   tImg : template image
%         tPos : template center
%         tDim : template size [nRows nCols]
%         sImg : search image
%         dmax : max displacement
%         subpix: method to etsimate subpixel displacement
%                 'none'
%                 'biquad' -> biquadradic fit around peak used
%         d0   : (optional) initial guess of displacement
% OUTPUT  d    : displacement estimate
%         ccMap: cross-correlation map
%         x1,x2: matrices with x1/x2-displacement values
%                useful to display ccMap as a mesh calling
%                mesh(x1,x2,ccMap)
% Liya with some range checkings

% set defaults
d = [0,0];
maxRow = d(1);
maxCol = d(2);
ccMap=[];
[nRowcCcMap,nColCcMap]=size(ccMap);
if(nargin<7)
   d0 = [0,0];
end;

% check whether tSize is odd
tOdd = mod(tDim,2);
if(~(tOdd(1) & tOdd(2)))
   error('template dimension is NOT odd');
end;

% crop the template
tWidth = (tDim(1)-1);
tHeight = (tDim(2)-1);
tRect = [tPos(1)-tWidth/2,tPos(2)-tHeight/2,tWidth,tHeight];
tPatch = imcropnint(tImg,tRect,'bilinear');
if(isempty(tPatch))
   error('template outside template image');
end;

% crop the search area
while( (maxRow < 1) | (maxRow > (nRowCcMap-1)) | ...
      (maxCol < 1) | (maxCol > (nColCcMap-1)))
	sPos = tPos + d0;
	sWidth = 2*ceil(abs(dmax(1)));
	sHeight = 2*ceil(abs(dmax(2)));
   if( (sWidth < 2 ) | (sHeight < 2) )
      error('dmax with non-zero components requested');
	end;
	sArHeight = sHeight + tHeight;
	sArWidth = sWidth + tWidth;
	sArRect = [sPos(1) - sArWidth/2,sPos(2) - sArHeight/2,...
   	   sArWidth,sArHeight];
   % Liya added here
   if sArRect(1)<=0
       sArRect(3) = sArRect(3) + sArRect(1) -1;
       sArRect(1) = 1;       
   end
   
   if sArRect(2)<=0
       sArRect(4) = sArRect(4) + sArRect(2) -1;
       sArRect(2) = 1;       
   end
   
   if sArRect(1)+ sArRect(3) > size(sImg,2)-1
      sArRect(3) = size(sImg,2)-2 - sArRect(1);
   end
      
   if sArRect(2)+ sArRect(4) > size(sImg,1)-1
      sArRect(4) = size(sImg,1)-2 - sArRect(2);
   end
   
   %
	sArPatch = imcropnint(sImg,sArRect,'bilinear');
	if( isempty(sArPatch))
   	error('search area outside search image');
	end;

	% compute cross-correlation
	ccMap = conv2(sArPatch,fliplr(flipud(tPatch)),'valid');

	% seek the maximum as an integer position
	[nRowCcMap,nColCcMap] = size(ccMap);
	rng = floor([nRowCcMap,nColCcMap]/2);
	[x1,x2] = meshgrid(-rng(2):rng(2),-rng(1):rng(1));
	[ccRowMax,ccRowMaxI] = max(ccMap);	
	[ccColMax,ccColMaxI] = max(ccRowMax);

	maxRow = ccRowMaxI(ccColMaxI);
	maxCol = ccColMaxI;
	d(1) = x1(maxRow,maxCol);
	d(2) = x2(maxRow,maxCol);

	% update d0
	% if the maximum is at the border of the currnet map shift
	% the maximum towards the center (is tested in the next while)
	d0 = d + d0;
end;

% determine the maximum in the cross-correlation map
if(strcmp(subpix,'none'))
   d = d0;
end;

if(strcmp(subpix,'biquad'))
   % generate 3x3 coordinate mask
   [x1,x2] = meshgrid(-1:1,-1:1);
   % crop 3x3 window around integer peak
   ccWin = ccMap(maxRow-1:maxRow+1,maxCol-1:maxCol+1);
   % fit around the peak
   [a,sa] = biquadfit(x1,x2,ccWin);
   denom = 4*a(4)*a(5) - a(6)^2;
   d(1) = (a(6)*a(3) - 2*a(2)*a(5));
   d(2) = (a(6)*a(2) - 2*a(3)*a(4));
   d = d/denom;
   d = d + d0;
end;

