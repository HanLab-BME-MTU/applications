function overlayTrack(r,p)

[fileName,dirName] = uigetfile('*.tif','Choose a .tif file');
I = imread([dirName,filesep,fileName]);
[xmax,ymax] = size(I);

if r == 5
    load([dirName(1:end-8),'\point_files\config001_5p00_track_bidir.mat']);
elseif r == 6
    load([dirName(1:end-8),'\point_files\config001_3p00_track_bidir.mat']);
end
indx = find( [tracks.len] >= 7);
le  =length(indx)
tracks = tracks(indx);
% find speeds
for i = 1 : le 
    d = diff(tracks(i).points(:, 1:2));
    tracks(i).meanDisp = sum(sqrt(sum(d.^2, 2)))/(tracks(i).len-1) ;
end
% find the fastest or the longest
if p == 1 % speed
    maxP = max([tracks.meanDisp]);
    indxP = find( [tracks.meanDisp]==maxP);
elseif p == 2 % length
    maxP = max([tracks.len]);
    indxP = find( [tracks.len]==maxP);
end
traj = tracks(indxP(1));

aaux = 5;
s = 3;
strg=sprintf('%%.%dd',s);
for i = 1:traj.len
    im = traj.startID + i -1; % get the current image nb
    indxStr=sprintf(strg,im);
    
    I = imread([dirName,fileName(1:end-7),indxStr,'.tif']);     
    If=Gauss2D(I,1);
    figure, imshow(If(1+aaux:end-aaux,1+aaux:end-aaux),[]);%I4
    hold on
    load([dirName(1:end-7),'cands\feats',indxStr],'feats')
    for j = 1:length(feats.ori)
        h = quiver(feats.pos(j,1)-aaux,feats.pos(j,2)-aaux,-cos(feats.ori(j)*pi/180),sin(feats.ori(j)*pi/180),3,'r');
        set(h,'LineWidth',2)
    end
    plot(traj.points(1:i,1)-aaux,traj.points(1:i,2)-aaux,'g-')
    plot(traj.points(i,1)-aaux,traj.points(i,2)-aaux,'g*')
end

traj

