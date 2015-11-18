function [ colors ] = getColors( data )
%getColors Get colors for each element in data

condColorAll = {...
    [0 0 0],... %black 1
    [0 0 1],... %blue 2
    [0 0.5 0],... %dark green 3
    [1 0 0],... %red 4
    [1 0.6 0.3],... %orange 5
    [0 1 0],... %green 6
    [1 0 1],... %magenta 7
    [0 1 1],... %cyan 8
    [1 1 0],... %yellow 9
    [0.7 0.5 0],... %brown 10
    [0.7 0.7 0.7]... %gray 11
    [0.5 0.5 1],... %purple 12
    [0.3 0.8 1],... %light blue 13
    };
colors = condColorAll(mod(1:numel(data), 13) + 1);
colors = reshape(colors,size(data));

end

