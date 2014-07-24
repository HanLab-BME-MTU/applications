function info=addTimeToFeatures(features,varargin)
%ADDTIMETOFEATURES adds time/frame information into last column of LM data
%
%  Input
%    features  ->  cell array as obtained from analyzeLMdata (required)
%    
%
%


ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('features',@iscell);

ip.parse(features,varargin{:});

f=ip.Results.features;

nf=numel(f);

info={};

for n=1:nf
    tmp=f{n};
    if( ~isempty(tmp) )
        rows=size(tmp,1);
        cols=size(tmp,2);
        sig=NaN(rows,cols+1);
        sig(:,1:end-1)=tmp;
        sig(:,end)=n;
    else
        sig=NaN;
    end
    info{n}=sig;
end

info=info';

