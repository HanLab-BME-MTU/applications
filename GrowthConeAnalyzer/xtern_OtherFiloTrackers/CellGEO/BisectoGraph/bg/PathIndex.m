function [PathInd,PathLng]=PathIndex(INDX,Vb)

% This function indexes each path from each boundary point

% ----------------- Input: ---------------------------


% ----------------- Output: ---------------------------


% -------------------------------------------------------

    Lb=length(Vb);
    path=zeros(1000,Lb);
    PathLng=zeros(1,Lb);
    for i=1:Lb    
        lnum=find(INDX(:,1)==Vb(i),1);
        p=0;
        while ~isempty(lnum)
            p=p+1;
            vnum=INDX(lnum,1);            
            path(p,i)=vnum;
            nnum=INDX(lnum,2);
            lnum=find(INDX(:,1)==nnum,1);
        end
        PathLng(i)=p+1;
        path(p+1,i)=INDX(1,2);
        
    end
    
    PathInd=path(1:max(PathLng),:);
    
    