function nAct=getUVActivation(features,varargin)
%GETUVACTIVATION determines number of signals due to UV activation pulse
%    
%  Input:
%    features  ->  cell array as obtained from analyzeLMdata (required)
%         nCh  ->  # of channels (optional, default 1)
%         rep  ->  # of frames in activation cycle (optional, default 10)
%


ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('features',@iscell);
ip.addOptional('nCh',1,@isscalar);
ip.addOptional('rep',10,@isscalar);

ip.parse(features,varargin{:});

F=ip.Results.features;
nCh=ip.Results.nCh;
rep=ip.Results.rep;

=numel(F);
nCyc=nf/(nCh*rep);

nAct=zeros(nCyc,nCh);

k=1;
for i=1:nCyc
    for c=1:nCh
        tmp=F{k};
        if( ~isempty(tmp) )
            nAct(i,c)=numel(tmp(:,1));
        else
            nAct(i,c)=0;
        end
        k=k+rep-1;
    end
end


end