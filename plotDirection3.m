function plotDirection3(linesIn,varargin)
%PLOTDIRECTION plots the input 3D lines(s) with an arrow indicating the direction of the curve
%
% plotDirections(linesIn)
% plotContours(linesIn,plotStyle1,plotStyle2,...)
%
% Description:
%   
%   Plots the input 3D lines ASSUMING THEY ARE IN IMAGE COORDINATES [y,x,z]
%   and shows the direction the points occur within the line by placing a
%   small arrow at the first point.
% 
% Input:
% 
%   linesIn - A single 3xM or Mx3 matrix containing the line to plot, or a
%             cell array of 2xM or Mx2 matrices if multiple lines are to be
%             plotted.
% 
%   plotStyleString - Optional. A string (or strings) specifying the
%                     style/color to use when plotting the lines (same as
%                     with the plot command)
%
% Hunter Elliott
% 2/2012
%

if ~iscell(linesIn)
    linesIn = {linesIn};
end

if nargin >  1
    plotArgs = varargin;
else
    plotArgs = {'b'};
end

hold on

%Convert all the lines to 3xM
needsTranspose = cellfun(@(x)(size(x,1) ~= 3),linesIn);
linesIn(needsTranspose) = cellfun(@(x)(x'),linesIn(needsTranspose),'UniformOutput',false);

%Plot the lines
cellfun(@(x)(plot3(x(2,:),x(1,:),x(3,:),plotArgs{:})),linesIn)

%Plot the first point with arrow
iNotPt = cellfun(@(x)(size(x,2)> 1),linesIn); %Make sure the line has 2 pts!
cellfun(@(x)(quiver3(x(2,1),x(1,1),x(3,1),x(2,2)-x(2,1),x(1,2)-x(1,1),x(3,2)-x(3,1),5,plotArgs{:},'MaxHeadSize',10)),linesIn(iNotPt));
cellfun(@(x)(plot3(x(2,1),x(1,1),x(3,1),plotArgs{:},'Marker','x')),linesIn(iNotPt));
