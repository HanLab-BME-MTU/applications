function [ux, uy] = elPostEval(fem,options,tout,xCoord,yCoord)
%ELASTICPOST  Post-process the solution and output the solution values on the
%             specified coordinates and time period.
%
% SYNOPSIS :
%    u = elasticPost(fem,options,tout,xCoord,yCoord)
%
% INPUT :
%    tout : A 1D array of time steps over which the solution is output.
%       For static problem, it can be set to [].
%    xCoord : X coordinates of the positions where we want the solution. It
%       can be either 1D or 2D array.
%    yCoord : Y coordinates of the positions where we want the solution. It
%       can be either 1D or 2D array.
%
%    Note: 'xCoord' and 'yCoord' must be either 1D or 2D at the same time. 
%       When they are both 1D arrays, a regular grid will be created.
%
% OUTPUT :
%    [ux, uy] : x and y component of the solution calculated for the time steps 
%       specified in 'tout' and the space coordinates, (xCoord, yCoord). The
%       size of 'ux or uy' is 'size(xCoord)-by-length(tout)'. 
%       When length(tout) == 1 or tout = [], the size of 'ux' or 'uy' is 
%       'size(xCoord)'.

% check 'xCoord' and 'yCoord':
if ~isempty(xCoord) & ~isempty(yCoord)
   if ~isnumeric(xCoord) | ~isnumeric(yCoord)
      error(['''xCoord'' and ''yCoord'' should be arrays of numerical ' ...
         'values.  See ELASTIC2D.']);
   elseif ndims(xCoord) ~= ndims(yCoord)
      error(['''xCoord'' and ''yCoord'' must be arrays of the same ' ...
         'dimensions.  See ELASTIC2D.']);
   elseif ndims(xCoord) > 2
      error(['''xCoord'' and ''yCoord'' are 1D or 2D arrays.  ' ...
         'See ELASTIC2D.']);
   else
      for k = 1:2
         if size(xCoord,k) ~= size(yCoord,k)
            error(['''xCoord'' and ''yCoord'' must have the same size ' ...
               'when they are 2D arrays.  See ELASTIC2D.']);
         end
      end
   end
elseif isempty(xCoord)+isempty(yCoord) == 1
   error(['''xCoord'' and ''yCoord'' must be empty arrays at the same ' ...
      'time.  See ELASTIC2D.']);
end

[ux,uy] = postinterp(fem,'ux','uy',xCoord,yCoord);
