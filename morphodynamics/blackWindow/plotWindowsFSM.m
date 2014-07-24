function plotWindowsFSM(axHandle,tag,winFileName,winColor)


% 
% plotWindowsFSM(axHandle,tag,winFileName,color)
% 
% Plots the window polygons contained in the .mat file specified by
% winFileName in the color winColor. The polygons use the unique tag tag
% and are plotted on the axes specified by axHandle.
% 
% Hunter Elliott, 6/2009
% 

if nargin < 4 || isempty(axHandle) || isempty(tag) || isempty(winColor)
    error('All inputs must be included!')
end

load(winFileName)

if ~exist('winPoly','var')
    error(strcat('Specified file "',winFileName,' does not contain window variable allWinPoly !! Check file!!!'))
end

%Determine the number of windows and bands.
[nBands,nWindows] = size(winPoly);


%loop through all windows
for m = 1:nBands
    for n = 1:nWindows                
        %Plot the current window if it exists
        if ~isempty(winPoly(m,n).outerBorder(:)) && ~isempty(winPoly(m,n).innerBorder(:))
            
                patch([winPoly(m,n).outerBorder(1,:) winPoly(m,n).innerBorder(1,end:-1:1) winPoly(m,n).outerBorder(1,1)],...
                 [winPoly(m,n).outerBorder(2,:) winPoly(m,n).innerBorder(2,end:-1:1) winPoly(m,n).outerBorder(2,1)],...
                 winColor,'EdgeColor',winColor,'FaceColor','none','Tag',tag,'Parent',axHandle)
            

        end    
    end
end
