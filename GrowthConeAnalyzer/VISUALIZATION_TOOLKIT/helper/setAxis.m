function [ h ] = setAxis(visible,percentScreen,fontSize) 
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(visible) 
    visible = 'off'; 
end

if nargin<2 || isempty(percentScreen)
    percentScreen = 0.75;
end 
if nargin<3 || isempty(fontSize) 
    fontSize = 20;
end 

h = fsFigure(percentScreen,'visible',visible);
 axes1 = axes('Parent',h,'FontWeight','bold','FontSize',fontSize, ...
     'FontName','Arial');
 hold(axes1,'all')
end

