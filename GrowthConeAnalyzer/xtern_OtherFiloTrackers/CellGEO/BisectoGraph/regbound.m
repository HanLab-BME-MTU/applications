function [fb,Status]=regbound(bnd)

    b=round(bnd);
    
    % --- Check the boundary matrix dimensions ----------------------------
    
    if size(b,2)>2 && size(b,1)==2
        b=b'; 
        Status=1;
    elseif size(b,2)==2 && size(b,1)>2
        Status=1;
    else
        fb=bnd;
        Status=0;
        return;
    end
    % ---------------------------------------------------------------------
  
    % --- Remove repeated boundary points ---------------------------------
    
    j=1; nb=b;
    for i=2:size(b,1)
        r=(b(i,1)-b(i-1,1))^2+(b(i,2)-b(i-1,2))^2;
        if r>1e-9
            j=j+1;
            nb(j,:)=b(i,:);
        end
    end
    b=nb(1:j,:);
    % ---------------------------------------------------------------------
    
    % --- Remove the end point if the same as the first one ---------------
    
    if b(1,1)==b(end,1) && b(1,2)==b(end,2) && size(b,1)>3
        b=b(1:(end-1),:);
        Status=1;       
    elseif (b(1,1)~=b(end,1) || b(1,2)~=b(end,2)) && size(b,1)>3
        Status=1;
    else
        fb=bnd;
        Status=0;
        return;
    end
    %----------------------------------------------------------------------

    % --- Check and fix the diraction of the boundary ---------------------
    y1=b(:,1);
    y2=[b(2:end,1); b(1,1)];
    x1=b(:,2);
    x2=[b(2:end,2); b(1,2)];
    A=sum((x2-x1).*(y2+y1));
    if A<0
        b=flipud(b);
    end
    % ---------------------------------------------------------------------    

    bs = b;
    bs(:,1)=b(:,1)-min(b(:,1))+6;
    bs(:,2)=b(:,2)-min(b(:,2))+6;
    xlim = max(bs(:,2))+5;
    ylim = max(bs(:,1))+5;

    nF = poly2mask(bs(:,2),bs(:,1),ylim,xlim);
    nF(ylim*(bs(:,2)-1)+bs(:,1))=1;

    % --- Remove one-pixel-thick appendages -------------------------------

    clF=nF;
    nXb=bs(:,2);
    nYb=bs(:,1);
    nur=1; N=length(nXb);

    while nur>0 && N>2    

        dXf=[nXb(2:N); nXb(1)]-nXb;
        dYf=[nYb(2:N); nYb(1)]-nYb;

        dXb=-[dXf(N); dXf(1:(N-1))];
        dYb=-[dYf(N); dYf(1:(N-1))];

        ind=find(abs(dXf-dXb)+abs(dYf-dYb)==0);
        nur=length(ind);    
        if nur>0        
            for i=1:nur
                clF(nYb(ind(i)),nXb(ind(i)))=0;
            end
            for i=1:nur
                nXb(ind(i))=Inf; nYb(ind(i))=Inf;
                if ind(i)==N
                    nXb(1)=Inf; nYb(1)=Inf;
                else
                    nXb(ind(i)+1)=Inf; nYb(ind(i)+1)=Inf;
                end
            end
            nXb=nXb(nXb~=Inf); 
            nYb=nYb(nYb~=Inf); 
        end
        N=length(nXb);

    end
    % ---------------------------------------------------------------------
    
    if N<3
        fb=bnd;
        Status=0;
        return;
    end
        
    % --- Dilate one-pixel-thick linker -----------------------------------

    fF=clF;
    ind=zeros(1,N); 
    j=0;
    for i=1:N
        d=abs(nXb(i)-nXb)+abs(nYb(i)-nYb);
        if sum(d==0)>1
            j=j+1;
            ind(j)=i;        
        end
    end
    for i=1:j
        if nXb(ind(i))==max(nXb) || nYb(ind(i))==max(nYb)
            fF(nYb(ind(i))-1,nXb(ind(i)))=1;
            fF(nYb(ind(i)),nXb(ind(i))-1)=1;
        else
            fF(nYb(ind(i))+1,nXb(ind(i)))=1;
            fF(nYb(ind(i)),nXb(ind(i))+1)=1;
        end
    end
    fbnd=bwboundaries(fF,'noholes');
    fb=fbnd{1};
    fb(:,1)=fb(:,1)+min(b(:,1))-6;
    fb(:,2)=fb(:,2)+min(b(:,2))-6;
    % ---------------------------------------------------------------------

    % --- Remove the end point if the same as the first one ---------------
    
    if fb(1,1)==fb(end,1) && fb(1,2)==fb(end,2) && size(fb,1)>3
        fb=fb(1:(end-1),:);
        Status=1;       
    elseif (fb(1,1)~=fb(end,1) || fb(1,2)~=fb(end,2)) && size(fb,1)>3
        Status=1;
    else
        fb=bnd;
        Status=0;
        return;
    end
    %----------------------------------------------------------------------