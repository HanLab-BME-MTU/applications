function coef = calFemDOFCoef(fem,XY,fv,sigma)
%calFemDOFCoef : Calculate the Degree Of Freedom coefficients for given
%                function values at sampling data points.
%
% SYNTAX : coef = calFemDOFCoef(fem,XY,fv)
%    Given the function values 'fv' at a set of sampling data points 'XY',
%    calculate the Degree Of Freedom coefficients in the fem structure 'fem'
%    so that 'fem' with the calculated 'ceof' can be used as the finite 
%    element reprsentation of the function.
%
%    For example, after we get 'coef', we can then evaluate the function at
%    any other points as 
%       val = postinterp(fem,'v',[Px;Py],'u',coef);
%
% INPUT :
%    fem   : The fem structure in FEMLAB. See help femstruct.
%    XY    : A 2-by-n matrix that specifies the coordinates of the sampling 
%            data points where the first row gives the x-coordinates and the 
%            second row gives the y-coordinates.
%    fv    : The function values at 'XY'. It is an m-by-n matrix where m is the
%            number of functions.
%    sigma : A regularization parameter. Pass [] for the default value 0. It
%            has to be nonegative.
%
% OUTPUT :
%    coef : An m-by-DOF matrix where m is the number of functions and DOF is
%           the Degree Of Freedom of the 'fem' structure.

if ~isnumeric(XY) | size(XY,1) ~= 2
   error(['The coordinates of the sampling data points (2nd argument) ' ...
      'is not correctly defined.']);
end

if ~isnumeric(fv)
   error('The function values (3rd argument) are not numeric.');
end

if size(fv,2) ~= size(XY,2)
   error(['The length of the function values do not match the number of ' ...
      'data points.']);
end

if isempty(sigma)
   sigma = 0;
end

if ~isnumeric(sigma) | sigma < 0
   error(['The regularization parameter (4th argument) has to be a ' ...
      'nonegative numerical value.']);
end

%exclude points that are outside of the mesh domain of 'fem'.
[is,pe] = postinterp(fem,XY);
XY(:,pe) = [];
fv(:,pe) = [];

numPts  = size(XY,2);
numFuns = size(fv,1);

%To calculate these coefficients, we first construct the matrix whose columns
% are the values of each finite element basis function at the sampling data
% points, 'XY'.

%Get the Degree Of Freedom.
DOF    = flngdof(fem); 
coefFS = zeros(DOF,1);

name = fem.dim;
A    = zeros(numPts,DOF);
for k = 1:DOF
   coefFS(k) = 1;
   bspF = postinterp(fem,name,XY,'u',coefFS);
   A(:,k) = bspF.';
   coefFS(k) = 0;
end

coef = zeros(DOF,numFuns);
for k = 1:numFuns
   coef(:,k) = (A.'*A+sigma*eye(DOF))\(A.'*fv(k,:).');
end

coef = coef.';
