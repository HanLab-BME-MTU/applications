function fx=gfun_vol_f_xcomp(x,y)
% Calculates the volume force on the cell at position x,y and returns it in
% units [N]?
% Input:        x,y might be matrices of the same size. 
% global input: traction force map.
% Output:       fx must have the same size as the input x,y and specify the
%               the x-component of the forces at position x,y.

% for testing:
%fx=x;

global globForce

globForce_xcomp=TriScatteredInterp(globForce.pos(:,1),globForce.pos(:,2),globForce.vec(:,1));
fx=globForce_xcomp(x,y);

% This important check takes less than 0.01sec:
badVal=isnan(fx);
if any(badVal)
    % set those values to 0:
    display(['!!!',num2str(sum(badVal)),' NaN values set to 0!!!'])
    fx(badVal)=0;
end