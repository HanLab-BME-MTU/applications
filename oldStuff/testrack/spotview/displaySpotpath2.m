function displaySpotpath(slist)

figure;
hold on;
grid on;
daspect([1 1 1]);

NSPOTS=4;

COL_LIST={'b', 'b', 'r','g'};
LINEWIDTH=[0.5  0.5 2 2];

[oc(1:NSPOTS).coor] = deal([]);
[oc(1:NSPOTS).tp] = deal(0);

for t= 1:length(slist)
    if ~isempty(slist(t).sp)
        spidxlist=[];
        [newC(1:NSPOTS).coor]=deal([]);
        
        for s=1:length(slist(t).sp)
            spidxlist(s) = getspidx(slist(t).sp(s).type);
            newC(spidxlist(s)).coor=slist(t).sp(s).cord;
        end;
        
        for s=1:NSPOTS
            color=COL_LIST{s};
            if (t-oc(s).tp)==1
                style='-';
            else
                style ='-.';
            end;
            
            if isempty(newC(s).coor)
                % fusing kin?
                if s==3 & ~isempty(newC(4).coor)
                    newC(s).coor=newC(4).coor;
                elseif s==4 & ~isempty(newC(3).coor)
                    newC(s).coor=newC(3).coor;
                end;
            end;
            
            if isempty(oc(s).coor)
                % splitting kin?
                if s==3 & ~isempty(oc(4).coor)
                    oc(s).coor=oc(4).coor;
                elseif s==4 & ~isempty(oc(3).coor)
                    oc(s).coor=oc(3).coor;
                else
                    oc(s).coor=newC(s).coor;
                end;
            end;
            
            
            
            line=[oc(s).coor;newC(s).coor];
            %         plot3(line(:,2),line(:,1),line(:,3),'Color',color,'LineStyle',style,'LineWidth',LINEWIDTH(s));
            %            arrow3(oc(s).coor,newC(s).coor,'Color',color,'LineStyle',style);
            
            
            if ~isempty(newC(s).coor) & (oc(s).coor~=newC(s).coor)
                fc=oc(s).coor.*[0.05 0.05 0.2]-[3 3 2];
                tc=newC(s).coor.*[0.05 0.05 0.2]-[3 3 2];
                arrow3(fc,tc,[color style]);
                oc(s).coor=newC(s).coor;
                oc(s).tp=t;
            end;
        end;
        k = waitforbuttonpress;
        key=get(gcf,'CurrentCharacter');
        if k & key=='c'
            delete(gca);
            hold on;
            rotate3d on;
            grid on;
            daspect([1 1 1]);
        end;
    end;
        if k & key=='x'
            break;
        end;
end
%light('Position',[60,60,60]), lighting gouraud; alpha(.8);
%light('Position',[-10 -10 -10],'Style','local');
rotate3d on;
daspect([1 1 1]);

function sidx=getspidx(type);
switch type
case 'spb1'
    sidx=1;
case {'spb','spb '}
    sidx=2;
case {'kin','kin '}
    sidx=3;
case 'kin1'
    sidx=4;
otherwise
    sidx=0;
end;