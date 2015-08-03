function [ h ] = setAxis(visible) 
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(visible) 
    visible = 'off'; 
end 

h = fsFigure(0.75,'visible',visible); 
 axes1 = axes('Parent',h,'FontWeight','bold','FontSize',16,...
     'FontName','Arial');
 hold(axes1,'all')
end

