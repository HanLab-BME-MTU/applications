function tracks=fillTrackGaps(tracks)
    if(~isempty(tracks))
    se=[zeros(1,tracks.numTimePoints()) 1 ones(1,tracks.numTimePoints())];
    for tIdx=1:length(tracks)
        gi=tracks(tIdx).gapMask;
        if(any(gi))
            copyIdx=1:tracks(tIdx).lifetime;
            copyIdx(gi)=0;
            copyIdx=imdilate(copyIdx,se);
            tracks(tIdx).x=tracks(tIdx).x(copyIdx);
            tracks(tIdx).y=tracks(tIdx).y(copyIdx);
            tracks(tIdx).z=tracks(tIdx).z(copyIdx);
        end
    end
    end