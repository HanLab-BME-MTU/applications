function cands2=fsmPrepCheckDistance(cands2,cands1)

% fsmPrepCheckDistance uses a look-up-table to verify if a new (secondary or tertiary) speckle are not closer than 
% theoretically possible to a primary (for tertiary also secondary) one; if a candidate is closer - it is discarded
%
% SYNOPSIS    cands2=fsmPrepCheckDistance(cands2,cands1)
%
% INPUT      cands2      :   cands for the new/secondary speckles
%            cands1      :   cands for the old/primary speckles
%
% OUTPUT     cands2      :   updated cands for the new/secondary speckles
%
%
% DEPENDENCES   fsmPrepCheckDistance uses { createDistanceMatrix }
%               fsmPrepCheckDistance is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002

indx1=find([cands1.status]==1);%new
indx2=find([cands2.status]==1);%new

% P2=zeros(length(cands2),1);
P2=zeros(length(indx2),1);
Q2=P2;

% aux=length(cands2);
aux=length(indx2);
for i=1:aux
%     P2(i)=cands2(i).Lmax(1);
%     Q2(i)=cands2(i).Lmax(2);
    P2(i)=cands2(indx2(i)).Lmax(1);
    Q2(i)=cands2(indx2(i)).Lmax(2);
end 

R2=cat(2,P2,Q2);

% P1=zeros(length(cands1),1);
P1=zeros(length(indx1),1);
Q1=P1;

% aux=length(cands1);
aux=length(indx1);
for i=1:aux
%     P1(i)=cands1(i).Lmax(1);
%     Q1(i)=cands1(i).Lmax(2);
    P1(i)=cands1(indx1(i)).Lmax(1);
    Q1(i)=cands1(indx1(i)).Lmax(2);
end 

R1=cat(2,P1,Q1);

% D=createDistanceMatrix(R2,R1);
% nula=1e-10;
D=createSparseDistanceMatrix(R2,R1,5.21);

RemovePosSec=[]; % initialization
count=0; % initialization
r=0; % initialization
 
for i=1:size(D,1)   
    
    rowI=D(i,:); % primary
%     minDistRow=find(rowI<=5.2); % Change BECAUSE of SPARSE MATRIXfind(rowI>nula); 
    minDistRow=find(rowI<=5.2 & rowI>0);
    if ~isempty(minDistRow)
        
%         if cands2(i).status==1 % secondary
            
%             aux=find([cands1(minDistRow).status]==1); % vector (primary)
%             minDistRow=minDistRow(aux);
%             if ~isempty(minDistRow)
                
%                 RelInt=[cands1(minDistRow).deltaI]./cands2(i).deltaI;
                RelInt=[cands1(indx1(minDistRow)).deltaI]./cands2(indx2(i)).deltaI;

                r=5.2-3.5.*exp(.5-RelInt); % MAX 5.2 % r(i)=(4.1814*RelInt)/(2.2159+RelInt)+2.5613; % lsqNonLin Fit % NOT USED HERE
                d=rowI(minDistRow);
%                 d=d+eps; % Prevent division-by-zero error for overlapping speckle positions
                res=r./d;
                
                if  any(res>=1) % i.e. r>d
                    count=count+1;
%                     RemovePosSec(count)=i;
                      RemovePosSec(count)=indx2(i);
                end
%             end
%         end
    end
end  

if ~isempty(RemovePosSec)
    cands2(RemovePosSec)=[];
end

 