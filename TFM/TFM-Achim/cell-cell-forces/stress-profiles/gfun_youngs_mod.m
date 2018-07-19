function E=gfun_youngs_mod(x,y)
% Calculates the Young's modulus at position x,y and returns it in units [Pa]
% Input:        x,y might be matrices of the same size. 
% global input: actin intensity information.
% Output: E must have the same size as the input x,y and specify the
% Young's modulus at position x,y.

% For testing:
% E=100;

global globYoung

% The first two arguments could be abolished, the FAST '*linear' option 
% only works if x and y are equally spaced!
E=interp2(globYoung.xmat,globYoung.ymat,globYoung.val,x,y,'*linear');
