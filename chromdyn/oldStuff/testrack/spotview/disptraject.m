function trac=disptraject(path)
%DISPRAJECT find and draw spot trajectory
%
% SYNOPSIS trac=spottraject(slist)
%
% INPUT path   : path struct from spottraject
%

% c: 14/09/01	dT

%CONST DEFINITIONS
load_sp_constants;

COLORS ={'r' [0 1 1] 'g' 'y' 'm' 'c'};

figure;
hold on;
rotate3d on;
grid on;
shift=[0 0 0];
origin=path(1).sp1(1,:);
v1=[1 0 0];
%first marker
plot3(path(1).sp1(1,2),path(1).sp1(1,1),path(1).sp1(1,3),'r+');
colgrad=1/length(path);
for t=2:length(path)
    shift=(path(t).sp1(1,:)-path(t-1).sp1(1,:));
    pause;
    mask=(path(t).sp1~=0);
    mshift=ones(size(path(t).sp1,1),1)*shift.*mask;
    path(t).sp1=path(t).sp1-mshift;
    v2=path(t).sp1(3,:)-origin;
    v2=v2/norm(v2);
    rotdir=cross(v2,v1);
    alpha=180*acos(dot(v1,v2))/pi;
    path(t).sp1=rot3vec(path(t).sp1,rotdir,alpha,origin).*mask;
    
    if t>2
        for i=1:min(size(path(t).sp1,1),size(path(t-1).sp1,1))
            if ((path(t-1).sp1(i,1)~=0) & (path(t).sp1(i,1)~=0))
                lcord=[path(t-1).sp1(i,:);path(t).sp1(i,:)];
                plot3(lcord(:,2),lcord(:,1),lcord(:,3),'Color',COLORS{i});
            end;
        end;
    end;
    COLORS{2}=max([0 0 0],COLORS{2}-colgrad);
end;
trac=path;
