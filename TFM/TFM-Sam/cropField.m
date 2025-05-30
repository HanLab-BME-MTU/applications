function [croppedField] = cropField(inputField, nRows, cardDir)
% CROPFIELD takes an input displacement or force field and crops it by the
% value nRows in the cardinal direction specified by cardDir
%
% Inputs:
%       inputField   = The displacement or force field to be cropped
%       nRows        = The number of rows to remove from the field
%       cardDir      = The side of the field to crop from (N, S, E, or W,
%           non case-sensitive, corresponding to the cardinal directions 
%           w.r.t the matrix)
%
% Outputs:
%       croppedField = The cropped displacement or force field
%
% Example:
%       [croppedField] = cropField(exampleField,5,'E')
%   This example would crop the exampleField by 5 rows from the right side 
%   of the matrix and produce the croppedField as an output.
%

arguments
    inputField {mustBeField}
    nRows {mustBePositive,mustBeInteger}
    cardDir {mustBeTextScalar}
end

paredX = unique(inputField.pos(:,1));
paredY = unique(inputField.pos(:,2));

if strcmp(cardDir,'N') || strcmp(cardDir,'n')
    ymax = paredY(end-(nRows-1));
    inputField.pos(inputField.pos(:,2) >= ymax,:) = [];
    inputField.vec(inputField.pos(:,2) >= ymax,:) = [];
elseif strcmp(cardDir,'S') || strcmp(cardDir,'s')
    ymin = paredY(nRows);
    inputField.pos(inputField.pos(:,2) <= ymin,:) = [];
    inputField.vec(inputField.pos(:,2) <= ymin,:) = [];
elseif strcmp(cardDir,'E') || strcmp(cardDir,'e')
    xmax = paredX(end-(nRows-1));
    inputField.pos(inputField.pos(:,1) >= xmax,:) = [];
    inputField.vec(inputField.pos(:,1) >= xmax,:) = [];
elseif strcmp(cardDir,'W') || strcmp(cardDir,'w')
    xmin = paredX(nRows);
    inputField.pos(inputField.pos(:,1) <= xmin,:) = [];
    inputField.vec(inputField.pos(:,1) <= xmin,:) = [];
else
    error('Invalid direction identifier, must be (N, S, E, or W) or (n, s, e, or w) string scalar')
end

croppedField = inputField;
