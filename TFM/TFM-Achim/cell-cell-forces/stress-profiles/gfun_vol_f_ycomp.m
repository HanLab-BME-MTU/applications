function fy=gfun_vol_f_ycomp(x,y)
% Calculates the y-component of the volume force on the cell at position x,y and returns it in
% units [N]?
% Input:        x,y might be matrices of the same size. 
% global input: traction force map.
% Output:       fy must have the same size as the input x,y and specify the
%               the y-component of the forces at position x,y.

% for testing:
global globForce

globForce_ycomp=TriScatteredInterp(globForce.pos(:,1),globForce.pos(:,2),globForce.vec(:,2));
fy=globForce_ycomp(x,y);

% This important check takes less than 0.01sec:
badVal=isnan(fy);
if any(badVal)
    % set those values to 0:
    % display(['!!!',num2str(sum(badVal)),' NaN values set to 0!!!'])
    fy(badVal)=0;
end