%[data] = selectConditionData(varargin) displays a GUI for selection of individual movies from a list
% If no input is provided, this function first loads data sets using loadConditionData.
%
% Input (optional)
%         data : movie data structure returned by loadConditionData()
%
% Output
%         data : selected subset of movie data structures

% Daniel Nunez, 2008 (last modified 06/05/2013 by F.A.)

function [data] = selectConditionData(data)

% Load data sets
if nargin<1 || isempty(data)
    data = loadConditionData();
    loadData = true;
    while loadData
        str = [];
        while ~any(strcmpi(str, {'y','n'}))
            str = input('Load additional ''condition'' folders? [y/n]: ', 's');
        end
        if strcmpi(str, 'y')
            [data] = [data loadConditionData()];
        else
            loadData = false;
        end
    end
end

idx = listSelectGUI({data.source}, [], 'move', 1);
data = data(idx);
