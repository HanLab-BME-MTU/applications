function [f,f2x] = extractF(features,varargin)
%EXTRACTF gets specified fields from cell array obtained by analyzeLMdata
%   The following field names are permissible:
%   ----------------------------------------------------
%       x        y           A       s       c
%       x_pstd   y_pstd      A_pstd  s_pstd  c_pstd
%       x_init   y_init
%       sigma_r  SE_sigma_r  RSS
%       pval_Ar  mask_Ar     hval_Ar  hval_AD  isPSF
%   ----------------------------------------------------
%
%   By default, the fields x, y, A, x_pstd, y_pstd, A_pstd and isPSF are
%   extracted and stored in appropriate variables.
%   Further fields should be provided by a cell array containing the
%   desired field names.
%
%   US, 2012/10/04
%

% fields to be extracted by default
f2x={'x','y','A','x_pstd','y_pstd','A_pstd','isPSF'};

ip=inputParser;
ip.CaseSensitive=true;
ip.StructExpand=true;

ip.addRequired('features',@iscell);
ip.addOptional('interval',[],@isnumeric);
ip.addOptional('fnames',[],@iscell);

ip.parse(features,varargin{:});

F=ip.Results.features;
interval=ip.Results.interval;

% determine interval
if isempty(interval)
    interval=[1,numel(F)];
else
    if numel(interval) == 1
        interval=[1,interval];
    end
end
        
if ~isempty( ip.Results.fnames )
    f2x=horzcat(f2x,ip.Results.fnames);
end

nFields=numel(f2x);

startFrame=interval(1);
endFrame=interval(2);

np=0;
for k=startFrame:endFrame
    if ~isempty( F{k} )
        np=np+numel(F{k}.x);
    end
end

f=NaN(np,numel(f2x)+1);

w1=1;
for k=startFrame:endFrame
    
    if ~isempty( F{k} )
        
        frame=ones(numel(F{k}.x),1)*k;
        w2=w1+numel(F{k}.x)-1;
        
        for kk=1:nFields
            tmp=F{k}.(f2x{kk});
            f(w1:w2,kk)=tmp';
        end
        f(w1:w2,kk+1)=frame;
        w1=w2+1;
    end
end
         
        

end

