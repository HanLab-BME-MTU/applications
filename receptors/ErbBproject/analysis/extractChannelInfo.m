function channel=extractChannelInfo(features,varargin)
%EXTRACTCHANNELINFO separates features from every channel
%
%  Input:
%    features  ->  cell array as obtained from analyzeLMdata (required)
%         nCh  ->  # of channels (optional, default 2)
%         rep  ->  # of frames in activation cycle (optional, default 10)
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('features',@iscell);
ip.addOptional('nCh',2,@isscalar);
ip.addOptional('rep',10,@isscalar);

ip.parse(features,varargin{:});

F=ip.Results.features;
nCh=ip.Results.nCh;
rep=ip.Results.rep;

% number of data sets per channel
npc=numel(F)/(nCh);

channel=repmat(cell(npc,1),1,nCh);

n=1;
for k=1:npc/10
    nn=(k-1)*rep;
    for kk=1:nCh
        for kkk=1:rep
            channel{nn+kkk,kk}=features{n};
            n=n+1;
        end
    end
end

end

