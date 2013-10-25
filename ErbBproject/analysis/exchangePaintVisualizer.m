function [handle]=exchangePaintVisualizer(PL, varargin)
% This function generates a merged exchange paint graph 

%default labels
labels = {'EGFR','ErbB2','ErbB3','IGF-1R','cMet'};

cmap = {'r','g','b','k','m','y','c'};


ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('PL',@iscell);

ip.addOptional('labels',labels, @iscell);
ip.addOptional('cmap',cmap, @iscell);

ip.parse(PL,varargin{:});

labels = ip.Results.labels;
cmap = ip.Results.cmap;

handle = figure;
hold;


for i=1:numel(PL)
    temp = PL{i}.com;
    scatter(temp(:,1),temp(:,2),[cmap{mod(i-1,7)+1},'.']);   
end

legend(labels);

end
