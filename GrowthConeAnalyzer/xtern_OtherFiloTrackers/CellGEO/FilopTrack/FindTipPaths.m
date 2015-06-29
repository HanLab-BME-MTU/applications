function [filo,pd,tuck]=FindTipPaths(PathInd,PathLng,Coor,DIST,Rad,Tip,CrRad,cnd)

    Lt=length(Tip);
    filo=cell(1,Lt);
    pd=zeros(1,Lt);
    tuck=zeros(1,Lt);

    for j=1:Lt
        tpath=PathInd(1:PathLng(Tip(j)),Tip(j));
        COOR=Coor(tpath,:);
        LENG=DIST(tpath);
        HIGH=Rad(tpath);
        
        ovun=(HIGH(2:end)-CrRad).*(HIGH(1:(end-1))-CrRad);
    
        if cnd == 0
            ibase=find(ovun<=0,1,'last');
        elseif cnd == 1
            ibase=find(ovun<=0,1,'first');
        elseif cnd == 2
            ii=find(ovun<=0);
            if length(ii)>2
                ibase=ii(3);
            elseif length(ii)==1
                ibase=ii;
            else
                ibase=[];
            end
        end
        
        if isempty(ibase)
           pd(j)=max(LENG);  
           filo{j}=COOR;
           tuck(j)=tpath(end-1);
        elseif HIGH(ibase)==CrRad
           pd(j)=LENG(ibase);
           filo{j}=COOR(1:ibase,:);
           tuck(j)=tpath(ibase);
        else            
           xbase=COOR(ibase,1)+(COOR(ibase+1,1)-COOR(ibase,1))*(CrRad-HIGH(ibase))/(HIGH(ibase+1)-HIGH(ibase)); 
           ybase=COOR(ibase,2)+(COOR(ibase+1,2)-COOR(ibase,2))*(CrRad-HIGH(ibase))/(HIGH(ibase+1)-HIGH(ibase)); 
           pd(j)=max(LENG)-LENG(ibase)+sqrt((xbase-COOR(ibase,1))^2+(ybase-COOR(ibase,2))^2);
           filo{j}=[COOR(1:ibase,:);[xbase ybase]];    
           tuck(j)=tpath(ibase);
        end
        
    end
    
