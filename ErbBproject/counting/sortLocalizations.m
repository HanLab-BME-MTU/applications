function binList=sortLocalizations(features,imSize,varargin)
%SORTLOCALIZATIONS sorts single-molecule localizations into discrete grid
%
%   required input arguments:
%       features -> cell array, localizations as output from analyzeLMdata
%         imSize -> image size in pixels (squared image assumed)
%
%   optional input arguments:
%        da -> lattice spacing in pixel, default 1 pixel
%
%   output:
%       binList -> particle list for each bin
%
%   US, 2012/11/20
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('features',@iscell);
ip.addRequired('imSize',@isscalar);
ip.addOptional('da',1,@isscalar);

ip.parse(features,imSize,varargin{:});

f=ip.Results.features;
imSize=ip.Results.imSize;
da=ip.Results.da;

binList=cell(imSize/da,imSize/da);

nFrames=numel(f);

for k=1:nFrames
    if ~isempty(f{k})
        nLoc=numel(f{k}.x);
        loc=[f{k}.x' f{k}.y'];
        
        for l=1:nLoc
            xk=ceil(loc(l,1)/da);
            yk=ceil(loc(l,2)/da);
            
            info=[f{k}.x(l) f{k}.y(l) f{k}.A(l) f{k}.x_pstd(l) f{k}.y_pstd(l) f{k}.A_pstd(l) k];
            content=binList{xk,yk};
            
            content=vertcat(content,info);
            
            binList{xk,yk}=content;
        end
    end
end

for k=1:imSize/da
    for l=1:imSize/da
        if size(binList{k,l},1) < 5                    
            binList{k,l} = [];
        end
    end
end

end