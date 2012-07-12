function tracks=trackSignals(features,varargin)
%TRACKSIGNALS connects signals in consecutive frames to improve accuracy
%   
%  Input:
%    features  ->  cell array as output from analyzeLMdata for a single 
%                  channel (required)
%         rep  ->  # of frames in activation cycle (optional, default: 10)
%
%  Output:
%    tracks: structure with fields
%           .startEnd  ->  start and end frame of each track
%          .shortPath  ->  tracks, w/o leading and trailing zeros
%           .normPath  ->  tracks, normalized to origin
%           .residual  ->  deviation from mean path
%            .nTracks  ->  number of tracks
%
%  Ulrich Schmidt, 05/08/2012

ip=InputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('features',@iscell);
ip.addOptional('rep',10,@isscalar);

ip.parse(features,varargin{:});

F=ip.Results.features;
rep=ip.Results.rep;

% determine maximum number of tracks
ppf=pointsPerFrame(F);
ntMax=sum(ppf(:));
% allocate memory for speed and convenience
prel=-ones(ntMax,numel(F));

% initialize point relations with first non-empty data set
pos1=F{1};
if ~isempty(pos1)
    id1=pos1(:,14) == 1 & pos1(:,16) == 0;
    pos1=pos1(id1,:);
end
k=0;
while isempty(pos1)
    k=k+1;
    pos1=F{k};
    if ~isempty(pos1)
        id1=pos1(:,14) == 1 & pos1(:,16) == 0;
        pos1=pos1(id1,:);
    end
end

if k == 0
    k=1;
end

id1=pos1(:,14) == 1 & pos1(:,16) == 0;
for n=1:sum(id1)
    prel(n,k)=n;
end

% initial number of tracks is given by number of initially detected signals
nTracks=sum(id1);

for n=1:numel(F)
    
    pos1=F{n};
    pos2=[];
    if n < numel(F)
        pos2=F{n+1};
    end
    
    % both coordinate vectors are empty, nothing to do
    if isempty(pos1) && isempty(pos2)
        continue;
    end
    % pos2 contains data -> beginning of new track(s) in frame n
    if isempty(pos1) && ~isempty(pos2)
        % how many new tracks?
        addT=size(pos2,1);
        for kk=1:addT
            prel(nTracks+kk,n+1)=kk;
        end
        nTracks=nTracks+addT;
        continue;
    end
    % pos1 and pos2 contain data -> continue existing track(s)
    if ~isempty(pos1) && ~isempty(pos2)
        np1=size(pos1,1);
        np2=size(pos2,1);
        % create distance matrix between points in pos1 and pos2
        d12=distMat2(pos1(:,1:2),pos2(:,1:2));
        % all distances larger critical radius are marked illegal
        d12( d12 > 1.0 )=-1;
        trackIndex=prel(1:nTracks,n);
        try
            [link12,link21]=lap(d12,-1,0,1);
            link12=link12(1:np1);
            link21=link21(1:np2);
        catch
            continue;
        end
        
        for kk=1:np1
            nextp=link12(kk);
            %tIndex=find( trackIndex == kk );
            if nextp <= np2
                tIndex=find( trackIndex == kk);
                prel(tIndex,n+1)=nextp;
            end
        end
        
        for kk=1:np2
            prevp=link21(kk);
            % no link exists -> begin new track
            if prevp > np1
                nTracks=nTracks+1;
                prel(nTracks,n+1)=kk;
            end
        end
    end             
end

% calculate length of each track
trackLength=zeros(nTracks,1);
for i=1:nTracks
    trackLength(i)=sum( prel(i,:) > -1 );
end

shortPath=repmat(struct('pos',[],'CoM',[],'length',[-1 -1 -1]),nTracks,1);
longPath=repmat(struct('pos',[],'length',[]),nTracks,1);

for k=1:nTracks
    for i=1:numel(F)
        j=prel(k,i);
        pos=F{i};
        if isempty(pos)
            continue;
        end
        if j ~= -1
            shortPath(k).pos=vertcat(shortPath(k).pos,pos(j,:));
            longPath(k).pos=vertcat(longPath(k).pos,pos(j,1:2));
            if shortPath(k).length(1) < 0
                shortPath(k).length(1)=i;
            end
            shortPath(k).length(2)=i;
        else
            longPath(k).pos=vertcat(longPath(k).pos,[NaN NaN]);
        end
    end
    shortPath(k).CoM=mean(shortPath(k).pos(:,1:2),1);
    shortPath(k).length(3)=size(shortPath(k).pos,1);
    longPath(k).length=shortPath(k).length;
end

% filter paths that are short than 2
id=false(nTracks,1);
for k=1:nTracks
    if shortPath(k).length(3) > 1
        id(k)=true;
    end
end

normPath=longPath(id);
for k=1:sum(id);
    start=normPath(k).length(1);
    len=normPath(k).length(3);
    % shift all positions by first position
    shift=normPath(k).pos(start,:);
    normPath(k).pos(start,:)=[NaN NaN];
    for kk=1:len-1
        normPath(k).pos(start+kk,:)=normPath(k).pos(start+kk,:)-shift;
    end    
end

allNormPaths=zeros(sum(id),numel(F),2);
for k=1:sum(id)
    allNormPaths(k,:,:)=normPath(k).pos;
end

shift=nanmean(allNormPaths);
shift=squeeze(shift);
    
% normPath=cell(nTracks,1);
% for k=1:nTracks
%     p=shortPath{k};
%     if ~isempty(p)
%         p1=p(1,1:2);
%         np=size(p,1);
%         for i=1:np
%             p2=p(i,1:2);
%             normPath{k}=vertcat(normPath{k},p2-p1);
%         end
%     end
% end

tracks.shortPath=shortPath;
tracks.longPath=longPath;
%tracks.startEnd=startEnd;
%tracks.normPath=normPath;
tracks.nTracks=nTracks;
tracks.prel=prel;
tracks.shift=shift;

end