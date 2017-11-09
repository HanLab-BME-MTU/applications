function [xbase,ybase,lbase,ibase]=FindBase(rpath,lwcrv,CrR,cnd)

    ovun=(lwcrv(2:end,2)-CrR).*(lwcrv(1:(end-1),2)-CrR);
    
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
       xbase=rpath(end,1);
       ybase=rpath(end,2);
       lbase=lwcrv(end,1);       
    elseif lwcrv(ibase,2)==CrR       
       xbase=rpath(ibase,1);
       ybase=rpath(ibase,2);  
       lbase=lwcrv(ibase,1);
    else
       
       xbase=rpath(ibase,1)+(rpath(ibase+1,1)-rpath(ibase,1))*(CrR-lwcrv(ibase,2))/(lwcrv(ibase+1,2)-lwcrv(ibase,2)); 
       ybase=rpath(ibase,2)+(rpath(ibase+1,2)-rpath(ibase,2))*(CrR-lwcrv(ibase,2))/(lwcrv(ibase+1,2)-lwcrv(ibase,2)); 
       lbase=lwcrv(ibase,1)+(lwcrv(ibase+1,1)-lwcrv(ibase,1))*(CrR-lwcrv(ibase,2))/(lwcrv(ibase+1,2)-lwcrv(ibase,2)); 
    end