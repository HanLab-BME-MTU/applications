function dispspot3D(stack,spots)
% DISPSPOT3D display one timeframe of 3D stack with projection of image

%test if already open
fh=findobj(0,'Type','Figure','Tag','DispSpot3D');
% create or bring to front
if isempty(fh)
    fh=figure('Name','Spots','Tag','DispSpot3D');
else
    figure(fh);
end;
%plot first timestep
cla;
hold on;
grid on;
axis ij;
for i=1:size(spots,2)
    if strcmp(spots(i).type,'spb')
        mark='-.';
    elseif strcmp(spots(i).type,'kin')
        mark='+:';
    else
        mark='r+';
    end;
    stem3(spots(i).cord(1),spots(i).cord(2),spots(i).cord(3),mark);
end;
rotate3d on;
psl=max(stack,[],3);
axis([0 size(psl,2) 0 size(psl,1)])
psl=(psl-mean(psl(:)));
psl=psl/max(psl(:));
[AZ,EL] = view;
warp(zeros(size(psl)),psl);
view(AZ,EL);




