function trac=spottraject(slist)
%SPOTTRAJECT find and draw spot trajectory
%
% SYNOPSIS spottraject(slist)
%
% INPUT img   : stack time series
%

% c: 14/09/01	dT

%CONST DEFINITIONS
load_sp_constants;
global PIXELSIZE_XY PIXELSIZE_Z;

COLORS ={'r' 'b' 'g' 'y' 'm' 'c'};
curcol=COLORS;
figure;
hold on;
rotate3d on;
grid on;
%pixel 2 micron
p2m=ones(10,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
%prep spots
nsp1=size(cat(1,slist(1).sp.cord),1);
path(1).sp1=cat(1,slist(1).sp.cord).*p2m(1:nsp1,:);

ct=2;
for t=2:length(slist)
    if~isempty(slist(ct).sp)
        nsp2=size(cat(1,slist(ct).sp.cord),1);
        sp2=cat(1,slist(ct).sp.cord).*p2m(1:nsp2,:);
        if(nsp2>=nsp1)
            dm = distMat2(path(ct-1).sp1(:,:),sp2);
            [sdm insdm]=sort(dm(:));
            [idm jdm]=ind2sub(size(dm),insdm);
            in2=[];
            %remove unsused spots (cord= 0 0 0)
            zin1=find(path(ct-1).sp1(:,1)==0);
            for i=1:length(zin1)
                idm=idm.*(idm~=zin1(i));
            end;
            %assign correspodings spots starting with shortest distance 
            for i=1:length(idm)
                if idm(i)~=0
                    if jdm(i)~=0;
                        in2(idm(i))=jdm(i);
                        idm=idm.*(idm~=idm(i));
                        jdm=jdm.*(jdm~=jdm(i));
                    end;
                end;
            end;
            %            in2=dsearchn(sp2,sp1);
            in1=find(path(ct-1).sp1(:,1));
            in2=nonzeros(in2);
        else
            dm = distMat2(sp2,path(ct-1).sp1(:,:));
            [sdm insdm]=sort(dm(:));
            [idm jdm]=ind2sub(size(dm),insdm);
            in1=[];
             %assign correspodings spots starting with shortest distance 
            for i=1:length(idm)
                if idm(i)~=0
                    if jdm(i)~=0;
                        in1(idm(i))=jdm(i);
                        idm=idm.*(idm~=idm(i));
                        jdm=jdm.*(jdm~=jdm(i));
                    end;
                end;
            end;
            in2=[1:nsp2];
        end;
        pause;
        for i=1:min(nsp1,nsp2)
            lcord=[path(ct-1).sp1(in1(i),:);sp2(in2(i),:)];
            if(~any(lcord==0))
                line(lcord(:,2),lcord(:,1),lcord(:,3),'Color',COLORS{in1(i)});
            end;
        end;
        path(ct).sp1(in1,:)=sp2(in2,:);
%         %set colors
%         lc=1:length(curcol);
%         for i=1:length(in1) ,lc=nonzeros(lc.*(lc~=in1(i)));end;
%         lc=[in1 lc'];
%         curcol=curcol(lc);
        % add new spots at end;
        if(nsp2>nsp1)
            for i=1:nsp2
                zin1=find(path(ct).sp1(:,1)==0);
                if ~any(in2==i)
                    if~isempty(zin1)
                        path(t).sp1(zin1(1),:)=sp2(i,:);
                    else
                        path(t).sp1=[path(ct).sp1(:,:);sp2(i,:)];
                    end
                    %curcol(end+1)=COLORS(length(curcol)+1);
                end;
            end;
        end;
        
        nsp1=length(nonzeros(path(ct).sp1(:,1)));
        ct=ct+1;
    end;
end;
trac=path;
