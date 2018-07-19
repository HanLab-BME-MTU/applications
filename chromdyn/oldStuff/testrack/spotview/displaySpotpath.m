function displaySpotpath(slist)

figure;
hold on;
rotate3d on;
grid on;

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
        color=COL_LIST{sidx};
            if (t-oc(sidx).tp)==1
                style='-';
            else
                style ='-.';
            end;
            
            if isempty(oc(sidx).coor)
                % splitting kin?
                if sidx==3 & ~isempty(oc(4).coor)
                    oc(sidx).coor=oc(4).coor;
                elseif sidx==4 & ~isempty(oc(3).coor)
                    oc(sidx).coor=oc(3).coor;
                else
                    oc(sidx).coor=slist(t).sp(s).cord;
                end;
            end;
            if t<length(slist) & sidx>2
                if isempty(slist(t+1).sp)
                
            line=[oc(sidx).coor;slist(t).sp(s).cord];
            plot3(line(:,2),line(:,1),line(:,3),'Color',color,'LineStyle',style,'LineWidth',LINEWIDTH(sidx));
            oc(sidx).coor=slist(t).sp(s).cord;
            oc(sidx).tp=t;
        end;
        k = waitforbuttonpress;
        key=get(gcf,'CurrentCharacter');
        if k & key=='c'
            delete(gca);
            hold on;
            rotate3d on;
            grid on;
        end;
        
    end;
end

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