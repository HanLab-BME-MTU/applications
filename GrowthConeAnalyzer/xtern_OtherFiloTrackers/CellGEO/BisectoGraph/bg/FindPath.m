function [path,filo,lbase]=FindPath(INDX,Coor,DIST,Rad,vnum,CrRad)

% This function indexes Voronoi graph and identifies protrusion tip.

% ----------------- Input: ---------------------------


% ----------------- Output: ---------------------------


% -------------------------------------------------------

    path=zeros(1000,5);
    lnum=find(INDX(:,1)==vnum,1);
    p=0;
    while ~isempty(lnum)
        p=p+1;
        vnum=INDX(lnum,1);        
        path(p,:)=[vnum Coor(vnum,1) Coor(vnum,2) DIST(vnum) Rad(vnum)];
        nnum=INDX(lnum,2);
        lnum=find(INDX(:,1)==nnum,1);

    end
    path=path(1:p,:);

    
    ovun=(path(2:end,5)-CrRad).*(path(1:(end-1),5)-CrRad);
    ibase=find(ovun<=0,1,'last');
    if isempty(ibase)
       lbase=max(path(:,4));  % added 11/01/12 
       filo=path;
    elseif path(ibase,5)==CrRad
       lbase=path(ibase,4);
       filo=path(1:ibase,:);
    else
       xbase=path(ibase,2)+(path(ibase+1,2)-path(ibase,2))*(CrRad-path(ibase,5))/(path(ibase+1,5)-path(ibase,5)); 
       ybase=path(ibase,3)+(path(ibase+1,3)-path(ibase,3))*(CrRad-path(ibase,5))/(path(ibase+1,5)-path(ibase,5)); 
       lbase=max(path(:,4))-path(ibase,4)+sqrt((xbase-path(ibase,2))^2+(ybase-path(ibase,3))^2);
       filo=[path(1:ibase,:);[0 xbase ybase (max(path(:,4))-lbase) CrRad]];           
    end

