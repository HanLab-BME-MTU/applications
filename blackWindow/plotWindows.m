function plotWindows(windowIn,stringIn,showNum)
%PLOTWINDOWS plots the input windows on the current axes
% 
% plotWindows(windowIn)
% plotWindows(windowIn,plotStyle)
%
% This function is for displaying the window geometry. It plots the window
% oultines and fills in their interior on the current figure axes. This can
% be used to overlay the windows on an image or to inspect their geometry.
% 
% Input 
%
%   windowIn - The cell-array containing the windows, as created with
%   getMaskWindows or getMovieWindows.m
% 
%   stringIn - A plot style option(s) to pass to the patch command,
%   determining the appearance of the windows when plotted. May also be a
%   cell-array specifying multiple options, or a cell-array-of-cell-arrays
%   the same size as the window array, specifying different plot styles for
%   each window. Optional. Default is {'r','FaceAlpha',.2}, which plots
%   windows filled in transparently with red, and with a black outline.
%
%   showNum - Integer. If greater than 0, the indices of each showNum-th
%   window will be plotted next to the window. That is, if showNum equals
%   5, then every 5th window will have it's location in the window cell
%   array plotted next to it. Optional. Default is 0 (no numbers). WARNING:
%   If you have lots of windows, enabling this option may drastically slow
%   down the plotting.
%
% Output:
%
%   The windows will be plotted as polygons on the current axes.
%
%Hunter Elliott
%Re-Written 5/2010
%

if nargin < 1 || isempty(windowIn)
    error('Come on, you have to at least input the window cell-array!')
end

if nargin < 2 || isempty(stringIn)
    stringIn = {'r','FaceAlpha',.2};
elseif ~iscell(stringIn)
    stringIn = {stringIn};
end

if ~iscell(stringIn{1})
    %We need to replicate this string for each window            
    stringIn = cellfun(@(x)(arrayfun(@(y)(stringIn),1:numel(x),'UniformOutput',false)),windowIn,'UniformOutput',false);    
end

if nargin < 3 || isempty(showNum)
    showNum = 0;
end

prevHold = ishold(gca);%Get hold state so we can restore it.
if ~prevHold
    %If hold wasn't on, clear the axis and then turn it on
    cla
    hold on
end

for j = 1:numel(windowIn)        
    for k = 1:numel(windowIn{j})        
        if ~isempty(windowIn{j}{k})                
            currWin = [windowIn{j}{k}{:}];
            if ~isempty(currWin)

                patch(currWin(1,:),currWin(2,:),stringIn{j}{k}{:});

                if showNum && mod(j,showNum)==0 && mod(k,showNum) == 0
                    text(currWin(1,1),currWin(2,1),[num2str(j) ',' num2str(k)])                       
                end                    
            end
        end
    end       
end    

if ~prevHold %Restore previous hold state
    hold off
end
axis image
axis ij
