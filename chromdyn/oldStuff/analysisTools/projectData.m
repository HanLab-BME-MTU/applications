function idlist=projectData(idlist)
%projects data from idlist for visualization and analysis purposes
%
%SYNOPSIS idlist=projectData(idlist)
%
%INPUT idlist
%
%OUTPUT idlist with new fields:
%           idlist(t).analysis.rotMat
%           idlist(t).analysis.spbDist
%           idlist(t).analysis.cenCenter
%           idlist(t).analysis.spbCentroid : pos. of spb1
%           idlist(t).analysis.cenCentroid
%           idlist(t).analysis.cenProj(1-2)
%           idlist(t).analysis.cenDist(1-2)
%           idlist(t).analysis.cenDev(1-2)
%           to be continued...
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input
if any (findstr([idlist(1).stats.labelcolor{1:end}],'?'))
    idlist=autoLabel(idlist);
end

%find number of colors
nColors=size(idlist(1).linklist,1);
nFrames=size(idlist,2);

switch nColors
    case 1
        disp('no spindle found')
        return
    case 2 %only 2 spb's
        %find spindle 
        %SPINDLE ENDS ARE CONSIDERED TO BE LABELED WITH IDLIST.STATS.LABELLIST(1)
        spbList=findSpindle(idlist);
        %build coordinate lists
        [spb1,spb2]=coordLists(idlist,spbList,nFrames);
        %calc spindle center
        spbCentroid=spb1;
        
        %generate reference system (v: vector, e: unit vector, n: norm)
        %spb vector
        v_spb=spb2-spb1;
        [n_spb,e_spb]=normList(v_spb);
        
        %e_2: perpendicular to e_1, [1,0,0]
        v_2=cross(e_spb,[ones(nFrames,1),zeros(nFrames,2)]);
        [n_2,e_2]=normList(v_2);
        
        %e_3 perpendicular on e_1 and e_2
        v_3=cross(e_spb,e_2,2);
        [n_3,e_3]=normList(v_3);
        
        %set up rotation matrix
        rotMat=mkRotMat(e_spb,e_2,e_3);
        %calculate cenCentroid
        cenCentroid=[];
        
        %save all data in idlist
        for t=1:nFrames
            idlist(t).analysis.rotMat=rotMat(:,:,t);
            idlist(t).analysis.spbDist=n_spb(t);
            idlist(t).analysis.spbCentroid=spbCentroid(t,:);
        end
        
    case 3
        %find spindle 
        %SPINDLE ENDS ARE CONSIDERED TO BE LABELED WITH IDLIST.STATS.LABELLIST(1)
        [spbList,cenList]=findSpindle(idlist);
        %build coordinate lists
        [spb1,spb2,cen1]=coordLists(idlist,[spbList;cenList],nFrames);
        %calc spindle center
        spbCentroid=spb1;
        %calc cenCenter
        cenCenter=cen1;
        
        %generate reference system (v: vector, e: unit vector, n: norm)
        %e_1: spb-vector
        v_spb=spb2-spb1;
        [n_spb,e_spb]=normList(v_spb);
        
        %e_2: perp(spb-vector,cenCenter-spbAxis)
        v_cenAx=perpVector(spb1,e_spb,cenCenter);
        %test if enough spots to do proper e_2, if not do same as spb-only
        %(alternatively: look for closest good e_2 t+1->-1->2->-2 etc)
        for t=1:nFrames
            if ~isempty(idlist(t).linklist)
                if length(idlist(t).spot)<3
                    v_cenAx(t,:)=cross(e_spb(t,:),[1,0,0],2);
                    cenCenter(t,:)=spbCentroid(t,:)-v_cenAx(t,:);
                end
            end
        end
        [n_cenAx,e_cenAx]=normList(v_cenAx);
        
        %e_3 perpendicular on e_1 and e_2
        v_3=cross(e_spb,e_cenAx,2);
        [n_3,e_3]=normList(v_3);
        
        %set up rotation matrix
        rotMat=mkRotMat(e_spb,e_cenAx,e_3);
        
        %calculate cenCentroid
        cenCentroid=cenCenter+v_cenAx;
        
        %calculate projection of cen position on spindle axis
        cenProj1=sum((cen1-spb1).*e_spb,2);
        
        %calculate distance between cen and axis
        cenDist1=n_cenAx;
        
        %calculate projection of cen position on e_3
        cenDev1=sum((cen1-spb1).*e_3,2);
        
        %save all data in idlist
        for t=1:nFrames
            idlist(t).analysis.rotMat=rotMat(:,:,t);
            idlist(t).analysis.spbDist=n_spb(t);
            idlist(t).analysis.cenCenter=cenCenter(t,:);
            idlist(t).analysis.spbCentroid=spbCentroid(t,:);
            idlist(t).analysis.cenCentroid=cenCentroid(t,:);
            idlist(t).analysis.cenProj=cenProj1(t);
            idlist(t).analysis.cenDist=cenDist1(t);
            idlist(t).analysis.cenDev=cenDev1(t);
        end
        
    case 4
        %find spindle 
        %SPINDLE ENDS ARE CONSIDERED TO BE LABELED WITH IDLIST.STATS.LABELLIST(1)
        [spbList,cenList]=findSpindle(idlist);
        %build coordinate lists
        [spb1,spb2,cen1,cen2]=coordLists(idlist,[spbList;cenList],nFrames);
        %calc spindle center
        spbCentroid=spb1;
        %calc cenCenter
        cenCenter=0.5*(cen1+cen2);
        
        %generate reference system (v: vector, e: unit vector, n: norm)
        v_spb=spb2-spb1;
        [n_spb,e_spb]=normList(v_spb);
        
        %e_2: perp(spb-vector,cenCenter-spbAxis)
        v_cenAx=perpVector(spb1,e_spb,cenCenter);
        %test if enough spots to do proper e_2, if not do same as spb-only
        %(alternatively: look for closest good e_2 t+1->-1->2->-2 etc)
        for t=1:nFrames
            if ~isempty(idlist(t).linklist)
                if length(idlist(t).spot)<3
                    v_cenAx(t,:)=cross(e_spb(t,:),[1,0,0],2)*0.000001;
                    cenCenter(t,:)=spbCentroid(t,:)-v_cenAx(t,:);
                end
            end
        end
        [n_cenAx,e_cenAx]=normList(v_cenAx);
        
        %%e_3 perpendicular on e_1 and e_2
        v_3=cross(e_spb,e_cenAx,2);
        [n_3,e_3]=normList(v_3);
        
        %set up rotation matrix
        rotMat=mkRotMat(e_spb,e_cenAx,e_3);
        
        %calculate cenCentroid
        cenCentroid=cenCenter+v_cenAx;
        
        %calculate projection of cen position on spindle axis
        cenProj1=sum((cen1-spb1).*e_spb,2);
        cenProj2=sum((cen2-spb1).*e_spb,2);
        
        %calculate distance between cen and axis
        cenV_1=perpVector(spb1,e_spb,cen1);
        cenDist1=normList(cenV_1);
        cenV_2=perpVector(spb1,e_spb,cen2);
        cenDist2=normList(cenV_2);
        
        %calculate projection of cen position on e_3
        cenDev1=sum((cen1-spb1).*e_3,2);
        cenDev2=sum((cen2-spb1).*e_3,2);
        
        %save all data in idlist
        for t=1:nFrames
            idlist(t).analysis.rotMat=rotMat(:,:,t);
            idlist(t).analysis.spbDist=n_spb(t);
            idlist(t).analysis.cenCenter=cenCenter(t,:);
            idlist(t).analysis.spbCentroid=spbCentroid(t,:);
            idlist(t).analysis.cenCentroid=cenCentroid(t,:);
            idlist(t).analysis.cenProj=[cenProj1(t),cenProj2(t)];
            idlist(t).analysis.cenDist=[cenDist1(t),cenDist2(t)];
            idlist(t).analysis.cenDev=[cenDev1(t),cenDev2(t)];
        end
end

%calc XCoordinates
cenNDR=zeros(nFrames,3);
for t=1:nFrames
    if ~isempty(idlist(t).linklist)
        cenNDR(t,:)=idlist(t).centroid;
    end
end
%no drift; rotated
idlist=calcXCoord(idlist,'XCoordNoDriftRot',cenNDR,rotMat);
%spb-corrected
idlist=calcXCoord(idlist,'XCoordSpb',spbCentroid);
%spb-corrected and rotated
idlist=calcXCoord(idlist,'XCoordSpbRot',spbCentroid,rotMat);


save(['pdat-',date]);