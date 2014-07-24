nTracks=numel(tracksFinal);

% get proper (p) and improper (ip) tracks
% proper tracks: start at beginning of activation cycle
% improper tracks: start in middle of activation cycle
tlenp=[];
tlenip=[];
idxp=[];
idxip=[];

for i=1:nTracks
    startFrame=tracksFinal(i).seqOfEvents(1);
    featIndx=tracksFinal(i).tracksFeatIndxCG;
    if( mod(startFrame,10) == 1 )
        tlenp=vertcat(tlenp,[numel(featIndx) startFrame]);
        idxp=vertcat(idxp,i);
    else
        tlenip=vertcat(tlenip,[numel(featIndx) startFrame]);
        idxip=vertcat(idxip,i);
    end
end

tracksp=tracksFinal(idxp);
allTracksp=cell(numel(tracksp),1);

for i=1:numel(tracksp)
    t=tracksp(i).tracksCoordAmpCG;
    [~,ncols]=size(t);
    t=reshape(t,8,ncols/8)';
    start=tracksp(i).seqOfEvents(1);
    stop=tracksp(i).seqOfEvents(2);
    t(:,end+1)=start:stop;
    allTracksp{i}=t;
end

tracksip=tracksFinal(idxip);
allTracksip=cell(numel(tracksip),1);
for i=1:numel(tracksip)
    t=tracksip(i).tracksCoordAmpCG;
    [~,ncols]=size(t);t=reshape(t,8,ncols/8)';
    start=tracksip(i).seqOfEvents(1);
    stop=tracksip(i).seqOfEvents(2);
    t(:,end+1)=start:stop;
    allTracksip{i}=t;
end

% transform cell arrays to matrices
posp=cell2mat(allTracksp);
posip=cell2mat(allTracksip);

% fit uncertainty vs. amplitude
[rp,mp,bp]=regression(log(posp(:,5))',log(posp(:,4))');
[rip,mip,bip]=regression(log(posip(:,5))',log(posip(:,4))');

% average track length vs. starting frame
nFrames=3000;
tstep=500;
aveLenp=zeros(nFrames/tstep,3);
i=1;
k=1;
t=k*tstep;
tmp=[];
while( 1 )
    tmp=horzcat(tmp,tlenp(i,1));
    % aveLenp(k,1)=aveLenp(k,1)+tlenp(i,1);
    % aveLenp(k,2)=aveLenp(k,2)+1;
    
    if( tlenp(i,2) > t )
        aveLenp(k,1)=t;
        aveLenp(k,2)=mean(tmp);
        aveLenp(k,3)=std(tmp);
        tmp=[];
        k=k+1;
        t=k*tstep;
    end
    
    i=i+1;
    if( i > size(tlenp,1) )
        aveLenp(k,1)=t;
        aveLenp(k,2)=mean(tmp);
        aveLenp(k,3)=std(tmp);
        break;
    end
    
    if( t > nFrames )
        break;
    end
    
    
end 
    
aveLenip=zeros(size(aveLenp));
i=1;
k=1;
t=k*tstep;
while( 1 )
    tmp=horzcat(tmp,tlenip(i,1));
    % aveLenip(k,1)=aveLenip(k,1)+tlenip(i,1);
    % aveLenip(k,2)=aveLenip(k,2)+1;
    
    
    if( tlenip(i,2) > t )
        aveLenip(k,1)=t;
        aveLenip(k,2)=mean(tmp);
        aveLenip(k,3)=std(tmp);
        tmp=[];
        k=k+1;
        t=k*tstep;
    end
    
    i=i+1;
    if( i > size(tlenip,1) )
        aveLenip(k,1)=t;
        aveLenip(k,2)=mean(tmp);
        aveLenip(k,3)=std(tmp);
        break;
    end
    
    if( t > nFrames )
        break;
    end
end 
    
    
    