function [INDX,Vb,Coor,DIST1,DIST2,Rad,Pr,smPr,Tip]=FindProfile(EdgeMatrix,Xb,Yb,Nconvex,Nconcave,Nedge,Smoo)

% This function indexes Voronoi graph and identifies protrusion tip.

% ----------------- Input: ---------------------------

% _Negde_ is the number of edges in the graph.

% _Nconvex_ is the number of convex boundary points.

% _Nconcave_ is the number of concave boundary points.
% Lb=Nconvex+Nconcave is the total number of boundary points.
% Nvert=2*(Nconvex+Nconcave)+1*Nconcave-2 is the number of graph vertices.

% _EdgeMatrix_ is an [Nedge x 4] matrix with each line [x1 y1 x2 y2]
% giving the coordinates of two vertices connected by an edge.

% _Xb_ and _Yb_ are Lb-long arrays of coordinates of the boundary 
% (simple polygonal chain), within which the graph was built.

% Smoo is a parameter for smoothing (by means of the Gaussian filter)
% the boundary profile. Smoo defines the width of the filter window.

% ----------------- Output: ---------------------------

% _INDX_ is an [Nedge x 2] index matrix, where the numbers [ch,ph] in each 
% line give indexes of the "child" and "parent" vertices of an edge.
% The graph-tree starts at the root, INDX(1,2), and ends on the boundary, 
% going from a parent vertex to a child one. Every parent has two children. 
% Every child has one parent accept for concave boundary points 
% that have two parents.

% _Vb_ is an Lb-long array of the boundary point indexes (graph vertex numbers).

% _Coor_ is [Nvert x 2] matrix which gives vertex coordinates.

% _DIST1_ and _DIST2_ are Nvert-long arrays of distances from 
% each vertex to the graph root along graph edges. 
% For all vertices accept concave boundary points DIST1=DIST2.

% _Rad_ is an Nvert-long array of minimal distances from each vertex
% to the boundary.

% _Pr_ is the boundary profile: it is an Lb-long array of distances 
% (path lengths) from each boundary point to the graph root along graph edges.
% For concave boundary points, the mean length of two possible paths is taken. 

% _smPr_ is the smoothed boundary profile (Lb-long).

% _Tip_ is an array of indexes of protrusion tips, so that 
% Xb(Tip) and Yb(Tip) are the coordinates of the tips.
% Each Tip corresponds to a maxima of Pr 
% on the interval between two local minima of the smPr.

% -------------------------------------------------------

    Nvert=2*(Nconvex+Nconcave)+1*Nconcave-2;
    Lb=length(Xb);
    RR=[EdgeMatrix(:,1:2); EdgeMatrix(:,3:4)];
    R=RR(:,1).^2+RR(:,2).^2;
    Rb=Xb.^2+Yb.^2;
    [U, Um, Un]=unique(R);
    Coor=RR(Um(1:Nvert),:);

    %disp('      Nvert length(U) Nedge Nconvex Nconcave Nconvex+Nconcave Lb');
    %disp([Nvert  length(U) Nedge Nconvex Nconcave Nconvex+Nconcave Lb]);

    Rad=zeros(1,Nvert);
    Cn=zeros(size(R));    % 1 if identified as parent, 0 otherwise
    Sn=zeros(1,Nvert);    % 2 if root, 1 if convex boundary point
                          % -1 if concave boundary point, 0 otherwise  
    %----------------------------------------------------------
    Mb=[max(Xb) min(Xb) max(Yb) min(Yb)];
    for i=1:Nvert
        Rad(i)=DistToBoundX(Coor(i,1),Coor(i,2),Xb,Yb,Mb);                
    end
    [vi vj]=max(Rad);
    Cand=vj; M=1; Sn(vj)=2;
    
    %----------------------------------------------------------
    INDX=zeros(Nedge,2);
    Vb=zeros(1,Lb);
    DIST1=zeros(Nvert,1);
    DIST2=zeros(Nvert,1);
    splt=zeros(Nvert,1);
    
    cnt=0; icnt=0;
    while M>0
        cnt=cnt+1;
        NewCand=zeros(size(Nvert));
        NewM=0;                
        for m=1:M
            vj=Cand(m);
            ph=find(vj==Un);
            if length(ph)==3
                Cn(vj)=1;
                ch=zeros(3,1);   
                for c=1:3
                    if ph(c)<=Nedge
                        ch(c)=ph(c)+Nedge;
                    else
                        ch(c)=ph(c)-Nedge;
                    end

                    if Cn(Un(ch(c)))==0
                        NewM=NewM+1;
                        icnt=icnt+1;
                        NewCand(NewM)=Un(ch(c));
                        INDX(icnt,1)=Un(ch(c));
                        INDX(icnt,2)=Un(ph(c));
                        DIST1(Un(ch(c)))=DIST1(Un(ph(c)))+sqrt((Coor(Un(ph(c)),1)-Coor(Un(ch(c)),1))^2+(Coor(Un(ph(c)),2)-Coor(Un(ch(c)),2))^2);
                    end                            
                end
            elseif length(ph)==2
                Sn(vj)=-1;
                Vb(find(U(vj)==Rb,1))=vj;                
            elseif length(ph)==1                        
                Sn(vj)= 1;
                Vb(find(U(vj)==Rb,1))=vj;
            else
                disp('Problem');
                return;                        
            end                                        
        end 
        M=NewM;
        Cand=NewCand(1:M);                
    end                
    DIST2=DIST1;
    
    tmp=find(Sn(INDX(1:icnt,1))==-1);
    %for i=1:Nedge
    %    if Sn(INDX(i,1))==-1
    for ti=1:length(tmp)    
        i=tmp(ti);
        if splt(INDX(i,1))==0
            splt(INDX(i,1))=1;
            DIST1(INDX(i,1))=DIST1(INDX(i,2))+sqrt((Coor(INDX(i,2),1)-Coor(INDX(i,1),1))^2+(Coor(INDX(i,2),2)-Coor(INDX(i,1),2))^2);
        else    
            DIST2(INDX(i,1))=DIST2(INDX(i,2))+sqrt((Coor(INDX(i,2),1)-Coor(INDX(i,1),1))^2+(Coor(INDX(i,2),2)-Coor(INDX(i,1),2))^2);
        end
    end
    %    end
    
    Pr1=DIST1(Vb); Pr2=DIST2(Vb);
    Pr=(Pr1+Pr2)/2;
    [Px,Py]=GFilterPer(1:Lb,Pr',Smoo);

    Pp=[Py(2:end), Py(1)];
    Pm=[Py(end), Py(1:(end-1))];
    Xdip=Px((Pp>Py) & (Pm>Py));
    Tip=zeros(size(Xdip));
    
    tt1=1:Xdip(1);
    tt2=Xdip(end):Lb;
    [vt1 it1]=max(Pr(tt1));
    [vt2 it2]=max(Pr(tt2));
    if vt1>=vt2
        Tip(1)=tt1(it1);
    else
        Tip(1)=tt2(it2);
    end
    for t=2:length(Xdip)
        tt=Xdip(t-1):Xdip(t);
        [vt it]=max(Pr(tt));        
        Tip(t)=tt(it);
    end

    smPr=Py;
