lambda=672;
kon=5;
koff=2;
kb=0.01;

x=16:16:240;
y=16:32:240;
pos=[];
for i=x
    for j=y
        pos=vertcat(pos,[i ,j]);
    end
end

% generate artifical STORM experiment
[stack,emState,settings]=generateSTORMmovie(lambda,kon,koff,kb,'pos',pos);

% analyze STORM experiment
f=analyzeSTORMmovie(stack,settings);
% get unique pixel values of localizations
nRep=settings.nRep;
xi=[];
yi=[];
for k=1:nRep:numel(f)
    xi=horzcat(xi,f{k}.x_init);
    yi=horzcat(yi,f{k}.y_init);
end

query=[pos(:,2) pos(:,1)];

% extract localizations
pest=extractF(f);
% perform KD ball query
idx=KDTreeBallQuery(pest(:,1:2),query,1);
% filter and sort localizations
plist=cell(size(idx));
for k=1:numel(idx)
    tmp=pest(idx{k},:);
    [~,sid]=sort(tmp(:,end));
    tmp=tmp(sid,:);
    plist{k}=tmp;
end

% get time trace of localizations in occupied bins
[tt,ft,e,erg]=getTimeTrace(plist,stack);