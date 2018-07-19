function [ h ] = setFigure( nx,ny,visible )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3 
    visible = 'off';
end 
%h = figure('Visible',visible); 
h = figure('Visible',visible,'Position',[50 50 nx ny]);
iptsetpref('ImshowBorder','tight');

set(h, 'InvertHardcopy', 'off');

set(h, 'PaperUnits', 'Points');

set(h, 'PaperSize', [nx ny]);

set(h, 'PaperPosition', [0 0 nx ny]); % very important
% 
%  set(h,'DefaultLineLineSmoothing','on');
%   set(h,'DefaultPatchLineSmoothing','on');

% Configure axes 
ha = axes('Position',[0 0 1 1],'Visible',visible); 


end

