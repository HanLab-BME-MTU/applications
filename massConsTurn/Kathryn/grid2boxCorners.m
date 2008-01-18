function [A]=grid2boxCorners(y,x)
%GRID2BOXCORNERS finds the area of quadrilaterals from an unstructured grid
%
% DESCRIPTION: grid2boxCorners uses two m x n matrices representing the x
% and y coordinates of points on a grid (regular or unstructured). It
% constructs a 5 x (m-1)*(n-1) matrix containing the vertices of
% quadrilaterals with vertices at the grid points. The 
%
% SYNOPSIS: [A]=grid2boxCorners(y,x)
%
% INPUT:  y : m x n matrix containing y-coordinates of grid (regular or not)
%         x : m x n matrix containing x-coordinates of grid (regular or not)
%
% OUTPUT: A : (m-1) x (n-1) matrix of the areas of the quadrilaterals
%             enclosed by the grid
%
% USERNAME: kathomps
%
% DATE: 20-Nov-2006


% example grids
% y=[1 1 1 ; 3 3 3 ; 4 4 4 ; 4.5 4.5 4.5]
% x=[1 3 4 ; 1 3 4 ; 1 3 4 ; 1 3 4.5 ]
% or
% y=[1 1 1 1.5 ; 3 3 3 3.5 ; 4 4 4 4.5 ]
% x=[1 3 4 4.5; 1 3 4 4.5; 1 3 4 4.5]
% imshow(zeros(500,500));
% hold on
% scatter(50*x(:),50*y(:),'r+')

% number of spaces (boxes) between n grid points is n-1
nBoxY=size(y,1)-1;
nBoxX=size(x,2)-1;

% all the top left coordinates from the grid
tlY=y(1:end-1,1:end-1);
tlX=x(1:end-1,1:end-1);

% all the bottom left coordinates from the grid
blY=y(2:end,1:end-1);
blX=x(2:end,1:end-1);

% all the bottom right coordinates from the grid
brY=y(2:end,2:end);
brX=x(2:end,2:end);

% all the top right coordinates from the grid
trY=y(1:end-1,2:end);
trX=x(1:end-1,2:end);

% x and y coordinates for box corners [tl; bl; br; tr; tl]
boxesX=[tlX(:)'; blX(:)'; brX(:)'; trX(:)'; tlX(:)'];
boxesY=[tlY(:)'; blY(:)'; brY(:)'; trY(:)'; tlY(:)'];

% find the area of each box
A=polyarea(boxesX,boxesY);

% test for case of self-intersection; give NaN for area
[bowtie]=testBowtie(boxesX, boxesY);
A(bowtie)=nan;

A=reshape(A,[nBoxY nBoxX]);