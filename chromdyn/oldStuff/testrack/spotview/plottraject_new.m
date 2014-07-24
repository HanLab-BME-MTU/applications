function plottraject_new(slist,idmap)
%PLOTTRAJECT plots path of spots global nearest neighbour


%CONST DEFINITIONS
load_sp_constants;
global PIXELSIZE_XY PIXELSIZE_Z;

%CONST DEFINITIONS
COLORS ={'r' 'b' 'g' 'y' 'm' 'c'};


%init vars /figures
curcol=COLORS;
colordef black;
figure;
hold on;
rotate3d on;
grid on;

%pixel 2 micron
p2m=ones(10,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];

for t=1:length(slist)-1
    if ~isempty(idmap(t).sp)
        ct=1;
        while((t+ct<length(idmap)) &isempty(idmap(t+ct).sp))
            ct=ct+1;
        end;
        
        nsp1=size(cat(1,slist(t).sp.cord),1);
        nsp2=size(cat(1,slist(t+ct).sp.cord),1);
        sp1=cat(1,slist(t).sp.cord).*p2m(1:nsp1,:);
        sp2=cat(1,slist(t+ct).sp.cord).*p2m(1:nsp2,:);
        
        k = waitforbuttonpress;
        key=get(gcf,'CurrentCharacter');
        if k & key=='c'
            delete(gca);
            grid on;
        end;
        
        for i=1:length(idmap(t).sp)
            toS=idmap(t).sp{i};
            for j=1:length(toS)
                lcord=[sp1(i,:);sp2(toS(j),:)];
                if(~any(lcord==0))
                    line(lcord(:,2),lcord(:,1),lcord(:,3),'Color',idmap(t).link(i));
                end;
            end;
        end;
        
        %swap colors
%                  tC=COLORS(1:length([idmap(t).sp{:}]));
%                  COLORS(1:length([idmap(t).sp{:}]))=COLORS([idmap(t).sp{:}]);
%                  COLORS([idmap(t).sp{:}])=tC;
    end;
end;
