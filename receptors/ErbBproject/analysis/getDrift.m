function [ drift, frameIndex] = getDrift(features,varargin)
%GETDRIFT determines stage drift of microscopy data
%
%
%   INPUT: 
%       features  -> cell array as obtained from analyzeLMdata
%   optional:
%       cutoff  -> maximum shift in pixels between two frames (default: 1)
%       amp     -> Sets a min amplitue for spots used in drift correction
%       (default = 0)
%   
%   OUTPUT:
%       drift  -> array with xy drift between frames
%       frameIndex  -> logical array indexing empty frames
%
% US, 2012/10/25
%
% Update 2012/12/31 JLW:
%        added amp feature
%        fixed a bug when there is only onetracking features in a frame
%         the code would set dx and dy equal to each other

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('features',@iscell);

ip.addOptional('cutoff',1.0,@isscalar);
ip.addOptional('amp',0.0,@isscalar);

ip.parse(features,varargin{:});

F=ip.Results.features;
rc=ip.Results.cutoff;
amp=ip.Results.amp;

nFrames=numel(F);

%applies amplitude criteria. If F.Amp <amp exclude points
if amp >0.0
    F = cellfun(@(x) ampSelector(x,amp),F);
end

drift=zeros(nFrames,2);
frameIndex=true(nFrames,1);

for i=1:nFrames-1
    
    if isempty( F{i} )
        frameIndex(i)=false;
        continue;
    else
        pos1=[F{i}.x' F{i}.y'];
        if isempty( F{i+1} )
            i=i+1;
        else
            pos2=[F{i+1}.x' F{i+1}.y'];
        end
    end
    
    % number of points in frame i and i+1
    np1=numel(pos1(:,1));
    np2=numel(pos2(:,1));
    
    % caluclate distance matrix and discard distances beyond cutoff
    dm=distMat2(pos1,pos2);
    dm(dm > rc )=-1;
    % link points between frame i and i+1
    try
        [l12,~]=lap(dm,-1,0,1);
    catch
        drift(i+1,1:2)=[0 0];
        continue;
    end
    
    l12=l12(1:np1);
    
    shift=[];
    % go through all points in set pos1[]
    for k=1:np1
        % check for matching point in set pos2[]
        j=l12(k);
        if j > np2
            continue;
        else
            p1=pos1(k,:);
            p2=pos2(j,:);
            % calculate difference
            shift=vertcat(shift,p2-p1);
        end
    end
    
    if isempty(shift)
        shift = NaN;
    end
    
    drift(i+1,1:2)=mean(shift,1);
end

end


%subbordinate function that processes the cell array to remove points with
%amplitudes larger than parameter amp

function [out]=ampSelector(c,amp)
    index=c.A > amp;
    out= {structfun(@(x)logicalSelector(x,index),c,'UniformOutput',false)};
end

function [out]=logicalSelector(S,index)
    S(index==0)=[];
    out = S;
end



