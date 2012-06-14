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
    id1=pos1(:,14) == 0 & pos1(:,16) == 1;
    pos1=pos1(id1,:);
end
k=1;
while isempty(pos1)
    k=k+1;
    pos1=F{k};
    if ~isempty(pos1)
        id1=pos1(:,14) == 0 & pos1(:,16) == 1;
        pos1=pos1(id1,:);
    end
end

id1=pos1(:,14) == 0 & pos1(:,16) == 1;
for n=1:sum(id1)
    prel(n,k)=n;
end

% initial number of tracks is given by number of initially detected signals
nTracks=0;

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
        d12( d12 > 100.0 )=-1;
        trackIndex=prel(1:nTracks,n);
        try
            [link12,link21]=lap(d12,-1,0,1);
            link12=link12(1:np1);
            link21=link21(1:np2);
            
            for kk=1:np1
                nextp=link12(kk);
                tIndex=find( trackIndex == kk );
                if nextp <= np2
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
        catch
            continue;
        end
    end             
end

% calculate length of each track
trackLength=zeros(nTracks,1);
for i=1:nTracks
    trackLength(i)=sum( prel(i,:) > -1 );
end

shortPath=repmat(struct('pos',[],'CoM',[],'length',[-1 -1 -1]),nTracks,1);

for k=1:nTracks
    for i=1:numel(F)
        j=prel(k,i);
        pos=F{i};
        if isempty(pos)
            continue;
        end
        if j ~= -1
            shortPath(k).pos=vertcat(shortPath(k).pos,pos(j,:));
            if shortPath(k).length(1) < 0
                shortPath(k).length(1)=i;
            end
            shortPath(k).length(2)=i;
        end
    end
    shortPath(k).CoM=mean(shortPath(k).pos(:,1:2));
    shortPath(k).length(3)=size(shortPath(k).pos,1);
end
            
% create normalized paths, i.e. start at origin
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
%tracks.startEnd=startEnd;
%tracks.normPath=normPath;
tracks.nTracks=nTracks;

end
% translate particle tracks to real world coordinates 
% shortPath=cell(nTracks,1);
% startEnd=-ones(nTracks,2);
% 
% for k=1:nTracks
%     for i=1:numel(F)
%         j=prel(k,i);
%         pos=F{i};
%         if isempty(pos)
%             continue;
%         end
%         %id=pos(:,14) == 0 & pos(:,16) == 1;
%         %pos=pos(id,:);
%         if j ~= -1
%             shortPath{k}=vertcat(shortPath{k},pos(j,:));
%             % set start & end frame
%             if startEnd(k,1) < 0
%                 startEnd(k,1)=i;
%             end
%             startEnd(k,2)=i;
%         end
%     end
% end

%    pos1=F{n};
%    pos2=[];
%    if n < numel(F)
%      pos2=F{n+1};
%    end
%    
%    if ~isempty(pos1) && isempty(pos2)
%        id1=pos1(:,14) == 0 & pos1(:,16) == 1;
%        pos1=pos1(id1,:);
%        if isempty(pos1)
%            continue;
%        else
%            pos3=[];
%            if n > 1
%                pos3=F{n-1};
%            end
%            % frame between two empty frames: tracks only one frame long
%            if isempty(pos3)
%                np1=size(pos1,1);
%                for kk=1:np1
%                    nTracks=nTracks+1;
%                    prel(nTracks,n)=kk;
%                end
%            end
%        end
%        continue;
%    end
%    
%    if isempty(pos1) && ~isempty(pos2)
%        id2=pos2(:,14) == 0 & pos2(:,16) == 1;
%        pos2=pos2(id2,:);
%        
%        if isempty(pos2)
%            continue;
%        end
%        np2=size(pos2,1);
%        for kk=1:np2
%            nTracks=nTracks+1;
%            prel(nTracks,n+1)=kk;
%        end
%        continue;
%    end
%    
%    if ~isempty(pos1) && ~isempty(pos2)
%        id1=pos1(:,14) == 0 & pos1(:,16) == 1;
%        pos1=pos1(id1,:);
%        id2=pos2(:,14) == 0 & pos2(:,16) == 1;
%        pos2=pos2(id2,:);
%        
%        if isempty(pos1) && isempty(pos2)
%            continue;
%        end
%        
%        np1=size(pos1,1);
%        np2=size(pos2,1);
%        % create distance matrix between points in pos1 and pos2
%        d12=distMat2(pos1(:,1:2),pos2(:,1:2));
%        % all distances larger than one pixel are marked illegal
%        d12( d12 > 100.0 )=-1;
%        try
%            [link12 link21]=lap(d12,-1,0,1);
%            
%            link12=link12(1:np1);
%            link21=link21(1:np2);
%            
%            for kk=1:np1
%                nextp=link12(kk);
%                tIndex=find( trackIndex == kk );
%                if nextp <= np2
%                    prel(tIndex,n+1)=nextp;
%                end
%            end
%            
%            for kk=1:np2
%                % get previous particle (if any)
%                prevp=link21(kk);
%                % new particle ? => create new path
%                if prevp > np1
%                    nTracks=nTracks+1;
%                    prel(nTracks,n+1)=kk;
%                end
%            end
%            trackIndex=prel(1:nTracks,n+1);
%        catch
%            continue;
%        end
%    end
%     pos1=F{n};
%     pos2=F{n+1};
%     
%     if isempty(pos1)
%         continue;
%     end
%     if isempty(pos2)
%         if isempty(pos3)
%             id1=pos1(:,14) == 0 & pos1(:,16) == 1;
%             pos1=pos1(id1,:);
%             if ~isempty(pos1);
%                 np1=size(pos1,1);
%                 for kk=1:np1
%                     nTracks=nTracks+1;
%                     prel(nTracks,n)=kk;
%                 end
%             end
%             pos3=pos1;
%         end
%         continue;
%     end
%     
%     id1=pos1(:,14) == 0 & pos1(:,16) == 1;
%     id2=pos2(:,14) == 0 & pos2(:,16) == 1;
%     pos1=pos1(id1,:);
%     pos2=pos2(id2,:);
%     
%     if isempty(pos1)
%         continue;
%     end
%     if isempty(pos2)
%         id1=pos1(:,14) == 0 & pos1(:,16) == 1;
%         pos1=pos1(id1,:);
%         if ~isempty(pos1);
%             np1=size(pos1,1);
%             for kk=1:np1
%                 nTracks=nTracks+1;
%                 prel(nTracks,n)=kk;
%             end
%         end
%         continue;
%     end
%     
%     if ~isempty(pos1) && ~isempty(pos2)
%         np1=size(pos1,1);
%         np2=size(pos2,1);
%         % create distance matrix between points in pos1 and pos2
%         d12=distMat2(pos1(:,1:2),pos2(:,1:2));
%         % all distances larger than one pixel are marked illegal
%         d12( d12 > 100.0 )=-1;
%         try
%             [link12 link21]=lap(d12,-1,0,1);
%             
%             link12=link12(1:np1);
%             link21=link21(1:np2);
%             
%             for kk=1:np1
%                 nextp=link12(kk);
%                 tIndex=find( trackIndex == kk );
%                 if nextp <= np2
%                     prel(tIndex,n+1)=nextp;
%                 end
%             end
%             
%             for kk=1:np2
%                 % get previous particle (if any)
%                 prevp=link21(kk);
%                 % new particle ? => create new path
%                 if prevp > np1
%                     nTracks=nTracks+1;
%                     prel(nTracks,n+1)=kk;
%                 end
%             end
%             
%             % update trackIndex
%             trackIndex=prel(1:nTracks,n+1);
%             
%         catch
%             continue;
%         end
%     else
%         continue;
%     end
%end