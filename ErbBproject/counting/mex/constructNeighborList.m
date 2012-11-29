%CONSTRUCTNEIGHBORLIST generates list of neighbors on a regular 2D grid
%
%   required input arguments:
%       grid -> grid size in x and y, square grid is assumed -> scalar
%        pbc -> periodic boundary condistions? true | false
%
%   output:
%       nbList -> cell array with neighbor list for each grid site
%
%   This function is implemented via a MEX-file.
%
%   US, 2012/11/21
%