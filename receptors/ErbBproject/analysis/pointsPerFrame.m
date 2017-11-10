function ppf=pointsPerFrame(features,varargin)
%POINTSPERFRAME gets number of detected signals per frame from LM features
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

N=numel(F)/(nCh*rep);

ppf=zeros(N,rep,nCh);

k=1;
for n=1:N
    for c=1:nCh
        for r=1:rep
            tmp=F{k};
            if( ~isempty(tmp) )
                id=tmp(:,14) == 0 & tmp(:,16) > 0;
                ppf(n,r,c)=numel(tmp(id,1));
            else
                ppf(n,r,c)=0;
            end            
            k=k+1;
        end
    end
end

end

