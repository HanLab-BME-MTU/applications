function [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,Coor,PathInd,PathLng,Neck)

    cnt=length(PeakInd);
    NeckInd=zeros(cnt,2);
    
    lp=cell(cnt,1);
    for j=1:cnt
        i=find(sum(PathInd==PeakInd(j))>0,1,'first');
        tpath=PathInd(1:PathLng(i),i);
        P=find(tpath==PeakInd(j),1,'first');
        
        lpath=tpath(P:end);                
        L=find(Rad(lpath)<Neck,1,'first');
        if ~isempty(L)
            NeckInd(j,1)=lpath(L); 
            NeckInd(j,2)=lpath(L-1);     
            lp{j}=lpath(L:end);
        end
    end
    lp=lp(NeckInd(:,1)>0);
    NeckInd=NeckInd(NeckInd(:,1)>0,:);
    
    [~,b,~]=unique(NeckInd(:,1));
    NeckInd=NeckInd(b,:);
    lp=lp(b);
    %NeckX=Coor(NeckInd(:,1),1)+(Coor(NeckInd(:,2),1)-Coor(NeckInd(:,1),1))*(Neck-Rad(NeckInd(:,1)))/(Rad(NeckInd(:,2))-Rad(NeckInd(:,1))); 
    %NeckY=Coor(NeckInd(:,1),2)+(Coor(NeckInd(:,2),2)-Coor(NeckInd(:,1),2))*(Neck-Rad(NeckInd(:,1)))/(Rad(NeckInd(:,2))-Rad(NeckInd(:,1))); 
    NeckNum=length(b);  
    
    ReckInd=zeros(NeckNum,2); 
    for i=1:NeckNum        
        
               
        D=Inf; nS=[];

        for s=1:cnt
            S=find(lp{i}==PeakInd(s),1,'first');
            if ~isempty(S)
                if S>1 && S-1<D
                    nS=S;
                    D=S-1;
                end                
            end
        end
               
        for s=1:NeckNum
            S=find(lp{i}==NeckInd(s,2),1,'first');
            if ~isempty(S)
                if S>1 && S-1<D
                    nS=S;
                    D=S-1;
                end                
            end
        end

        if ~isempty(nS)        
            spath=lp{i}(1:nS);                
            L=find(Rad(spath)<Neck,1,'last');
            if ~isempty(L)                
                ReckInd(i,1)=spath(L+1); 
                ReckInd(i,2)=spath(L);         
            end
        end

    end

    
    ReckInd=ReckInd(ReckInd(:,1)>0,:);
    [~,b,~]=unique(ReckInd(:,1));
    ReckInd=ReckInd(b,:);
    %ReckX=Coor(ReckInd(:,1),1)+(Coor(ReckInd(:,2),1)-Coor(ReckInd(:,1),1))*(Neck-Rad(ReckInd(:,1)))/(Rad(ReckInd(:,2))-Rad(ReckInd(:,1))); 
    %ReckY=Coor(ReckInd(:,1),2)+(Coor(ReckInd(:,2),2)-Coor(ReckInd(:,1),2))*(Neck-Rad(ReckInd(:,1)))/(Rad(ReckInd(:,2))-Rad(ReckInd(:,1))); 
    ReckNum=length(b);   

    ConeBnd=cell(1,NeckNum+1);
    A=zeros(1,NeckNum+1);
    cx=zeros(1,NeckNum+1);
    cy=zeros(1,NeckNum+1);
    
    BNDind=PathInd(1,:);
    tPathInd=PathInd;
    for i=1:ReckNum        
        t=find(sum(tPathInd==ReckInd(i,1))==0);
        tPathInd=tPathInd(:,t);
        BNDind=BNDind(t);
    end
    x=Coor(BNDind,1);
    y=Coor(BNDind,2);
    ConeBnd{1}=[x y];        
    %A(1)=polyarea(x,y);
    %cx(1)=mean(x);
    %cy(1)=mean(y);   
    [cx(1),cy(1),A(1)]=myCentroid(x,y);
    
    for j=1:NeckNum
        t=find(sum(PathInd==NeckInd(j,2))>0);
        x=zeros(length(t),1);
        y=zeros(length(t),1);
        tcnt=0;
        for i=1:length(t)
            tpath=PathInd(1:PathLng(t(i)),t(i));
            P=find(tpath==NeckInd(j,2), 1, 'first');
            hpath=tpath(1:P);
            ok=1;
            for r=1:ReckNum
                if sum(hpath==ReckInd(r,2))>0
                    ok=0;
                end
            end
            
            if ok        
                tcnt=tcnt+1;
                x(tcnt)=Coor(hpath(1),1);
                y(tcnt)=Coor(hpath(1),2);                
            end
        end
        ConeBnd{j+1}=[x(1:tcnt) y(1:tcnt)];        
        %A(j+1)=polyarea(x(1:tcnt),y(1:tcnt));
        %cx(j+1)=mean(x(1:tcnt));
        %cy(j+1)=mean(y(1:tcnt));        
        [cx(j+1),cy(j+1),A(j+1)]=myCentroid(x(1:tcnt),y(1:tcnt));
    end