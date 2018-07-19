function [spd,varargout] = spFlowAvgSpeed(flowField,dirV)
%spFlowAvgSpeed : Calculate the average speed of a flow field along one given
%                 direction or the normal a given curve.
%
% SYNTAX : 
%    spd = spFlowAvgSpeed(flowField,dirV)
%    [spd,pp] = spFlowAvgSpeed(flowField,dirV)
%
% INPUT :
%    flowField : An m-by-4 matrix where m is the number of flow vectors. Each
%                row gives the flow vector in the form [y0,x0,y1,x1] where
%                (y0,x0) is the base of the vector and (y1,x1) is the tip of
%                the vector.
%    dirV      : An n-by-2 matrix. Each row of the matrix has the form (y,x).
%                When n==1, it is the y- and x- coordinates of the given
%                direction. When n==2, it represents a line connected by the
%                two points given in the two rows and the speed is calculated
%                along the normal direction to the line. When n>2, a spline
%                curve is fitted throught the n points and the speed is
%                calculated along the normal directions to the curve. For each
%                flow vector, we first find the closest point on the curve and
%                calculate the speed of the flow vector along the normal at 
%                that point.
%
% OUTPUT :
%    spd : The average speed.
%    pp  : The pp-form of the cubic spline can be returned when 'dirV' is
%          list of points (>2).

%Check input.
if nargout > 2
   error('Too many input arguments.');
end

if ~isnumeric(flowField) | size(flowField,2) ~= 4
   error(['The first input is not in correct format. ' ...
      'See help spFlowAvgSpeed.']);
end

if ~isnumeric(dirV) | size(dirV,2) ~= 2
   error(['The second input is not in correct format. ' ...
      'See help spFlowAvgSpeed.']);
end

%Get the flow vector
numFlowVectors = size(flowField,1);
flowV          = flowField(:,3:4)-flowField(:,1:2);

%Check and parse the second input 'dirV'.
if size(dirV,1) == 1
   if nargout > 1
      error('Only one output when size(dirV,1)==1.');
   end
   dirV = dirV.'/norm(dirV);
   spd  = sum(flowV*dirV)/numFlowVectors;
elseif size(dirV,1) == 2
   if nargout > 1
      error('Only one output when size(dirV,1)==2.');
   end
   dirV = dirV(2,:)-dirV(1,:);
   dirV = [-dirV(2);dirV(1)]/norm(dirV);
   spd  = sum(flowV*dirV)/numFlowVectors;
else
   curvLen = diff(dirV);
   curvLen = sum(sqrt(curvLen.*curvLen),2);
   curvLen = [0;cumsum(curvLen)];
   curvLen = curvLen/curvLen(end);

   %Interpolate points in 'dirV' by cubic spline using 'curvLen' as parameter.
   pp  = csape(curvLen,dirV.');
   dpp = fnder(pp);

   %For each flow vector, identify the closest point on the curve.
   normal = zeros(numFlowVectors,2);
   for k = 1:numFlowVectors
      %First calculate the distance to each point of 'dirV'.
      distY = flowField(k,1)-dirV(:,1);
      distX = flowField(k,2)-dirV(:,2);
      dist = sqrt(distY.^2+distX.^2);
      [minD,ind] = min(dist);

      %We search for the closest point in the interval [leftS,rightS];
      if ind == 1
         leftS  = curvLen(1);
         rightS = curvLen(2);
      elseif ind == length(curvLen)
         leftS  = curvLen(end-1);
         rightS = curvLen(end);
      else
         leftS  = curvLen(ind-1);
         rightS = curvLen(ind+1);
      end

      [minS,minD] = fminbnd(@distFun,leftS,rightS,[],pp,flowField(k,1:2));
      normal(k,:) = fnval(dpp,minS).';
      normal(k,:) = [-normal(k,2),normal(k,1)]/norm(normal(k,:));
   end
   spd = sum(sum(flowV.*normal,2))/numFlowVectors;

   if nargout == 2
      varargout{1} = pp;
   end
end

function dist = distFun(s,pp,pt)
%distFun : Calculate the distance from 'pt' to the point on the curve given by
% parameter s and cubic spline interpolated curve.

curvPt = fnval(pp,s).';

dist = norm(pt-curvPt);
