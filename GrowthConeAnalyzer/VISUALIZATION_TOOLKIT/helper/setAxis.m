function [ h ] = setAxis
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

h = fsFigure(0.75,'visible','on'); 
 axes1 = axes('Parent',h,'FontWeight','bold','FontSize',16,...
     'FontName','Arial');
 hold(axes1,'all')
end

