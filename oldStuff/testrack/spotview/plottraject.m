function plottraject(slist)
%PLOTTRAJECT plots path of spots global nearest neighbour


%CONST DEFINITIONS
load_sp_constants;
global PIXELSIZE_XY PIXELSIZE_Z;
COLORS ={'r' 'b' 'g' 'y' 'm' 'c'};
%init vars /figures
curcol=COLORS;
colordef black;
ct=1;       %counter
figure;
hold on;
rotate3d on;
grid on;
%pixel 2 micron
p2m=ones(10,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
%prep spots
nsp1=size(cat(1,slist(1).sp.cord),1);
path(1).sp=cat(1,slist(1).sp.cord).*p2m(1:nsp1,:);
ct=ct+1;
while(nsp1>4)
    nsp1=size(cat(1,slist(ct).sp.cord),1);
    path(1).sp=cat(1,slist(ct).sp.cord).*p2m(1:nsp1,:);
    ct=ct+1;
end;

for t=2:length(slist)
    if ~isempty(slist(t).sp)
        nsp2=size(cat(1,slist(t).sp.cord),1);
        %make a switch statement with combinations of two vars (nsp1,nsp2)
        nsp=10*nsp1+nsp2;
        switch nsp
        case {22,33,34,43,44}     %(3 ->3, 3->4, 4 ->3, 4->4)
            sp2=cat(1,slist(t).sp.cord).*p2m(1:nsp2,:);
            [indSp1, indSp2, dist]=mappoints(path(ct-1).sp,sp2);
            
        case {42,32}   %(4 ->2, 3->2)
            sp2=cat(1,slist(t).sp.cord).*p2m(1:nsp2,:);
            % map only two closest         
            [indSp1, indSp2, dist]=mappoints(path(ct-1).sp,sp2);
            %             if isempty(indSp1)
            %                 break;
            %             end;
            [dist indist]=sort(dist);
            indSp1=indSp1(indist(1:2));
            indSp2=indSp2(indist(1:2));
            dist=dist(1:2);
            
        case {23,24}   %(2->3, 2->4)
            sp2=cat(1,slist(t).sp.cord).*p2m(1:nsp2,:);
            [indSp1, indSp2, dist]=mappoints(path(ct-1).sp,sp2);
            [dist indist]=sort(dist);

            path(ct-1).sp(3:length(indSp2),:)=sp2(indSp2(indist(3:end)),:);
            [indSp1, indSp2, dist]=mappoints(path(ct-1).sp,sp2);
%             indSp1=indSp1(indist(1:2));
%             indSp2=indSp2(indist(1:2));
        otherwise
            indSp1=[];
        end;
        %check for possible miss identi
        if ~isempty(indSp1)
            
        tl=find((dist<0.5 | dist<2*mean(dist(:))));
        indSp1=indSp1(tl);
        indSp2=indSp2(tl);
        
        path(ct).sp=sp2(indSp2,:);        
        %pause;
        for i=1:length(indSp1)
            lcord=[path(ct-1).sp(indSp1(i),:);sp2(indSp2(i),:)];
            if(~any(lcord==0))
                line(lcord(:,2),lcord(:,1),lcord(:,3),'Color',COLORS{i});
            end;
        end;
        nsp1=size(path(ct).sp,1);
        ct=ct+1;
    end;
        %         if(nsp2<5 & abs(nsp2-nsp1)<=1 &nsp2>1)
        %             sp2=cat(1,slist(t).sp.cord).*p2m(1:nsp2,:);
        %             if (nsp2<=nsp1)
        %                 indN=mappoints(path(ct-1).sp,sp2);
        %                 path(ct).sp=sp2(indN,:);
        %                 %path(ct).color=path(ct-1).color{};
        %                 pause;
        %                 for i=1:nsp1
        %                     lcord=[path(ct-1).sp(i,:);sp2(indN(i),:)];
        %                     if(~any(lcord==0))
        %                         line(lcord(:,2),lcord(:,1),lcord(:,3),'Color',COLORS{i});
        %                     end;
        %                 end;
        %             else
        %                 indN=mappoints(sp2,path(ct-1).sp);
        %                 [s invN]=sort(indN);
        %                 path(ct).sp=sp2(invN,:);
        %                 pause;
        %                 for i=1:nsp2
        %                     lcord=[path(ct-1).sp(s(i),:);sp2(invN(i),:)];
        %                     if(~any(lcord==0))
        %                         line(lcord(:,2),lcord(:,1),lcord(:,3),'Color',COLORS{i});
        %                     end;
        %                 end;
        %             end;
        %         end;
    end;
end;