function selectTracks(axis,traj,I,n)

if nargin == 0
    [fileName,dirName] = uigetfile('*.tif','Choose an image');
    I = imread([dirName,fileName]);
    load([dirName(1:end-8),'\point_files\config001_6p00_track_bidir.mat']);
    load([dirName(1:end-8),'\poles\axis.mat']);
    traj=tracks(find([tracks.len]>=7));
    n = 1;
end
leTraj = length(traj);
% calculate poles information for the n-th set of pole coordinates
dYp = poleAxis(n,3) - poleAxis(n,1);
dXp = poleAxis(n,4) - poleAxis(n,2);
poleVec = [dXp; dYp];
for i = 1:leTraj
    dY = traj(i).points(end,1) - traj(i).points(1,1);
    dX = traj(i).points(end,2) - traj(i).points(1,2);
    traj(i).vec = [dX; dY];
    traj(i).vel = sqrt(dY^2+dX^2)/traj(i).len;
    % cluster tracks based on distance to each pole of track's beginning
    dy1 = traj(i).points(1,1) - poleAxis(n,1);
    dx1 = traj(i).points(1,2) - poleAxis(n,2);
    dy2 = traj(i).points(1,1) - poleAxis(n,3);
    dx2 = traj(i).points(1,2) - poleAxis(n,4);
    dist1 = sqrt(dy1^2+dx1^2);
    dist2 = sqrt(dy2^2+dx2^2);
    if dist1 < dist2
        traj(i).pol = 2;
    else
        traj(i).pol = 1;
    end
end
% exclude tracks with NO motion
traj = traj(find([traj.vel]>0)); % SOME DONT MOVE!?
% traj = traj(find([traj.pol]==2)); % look at pole2 
traj = traj(find([traj.pol]==1)); % look at pole1
% calculate the dot product
magPoleVec = sqrt(sum(poleVec.^2,1));
magTrajVec = sqrt(sum([traj.vec]'.^2,2));
cos_trackAng = (poleVec'*[traj.vec])./(magPoleVec*magTrajVec)';
% select the positive values
% list = find(cos_trackAng<0); % for pole2
list = find(cos_trackAng>0); % for pole1
selTr = traj(list);
le = length(list);
% plot the poles and the tracks
figure,imshow(I,[])
hold on
plot(poleAxis(n,1),poleAxis(n,2),'b*')
plot(poleAxis(n,3),poleAxis(n,4),'g*')
for i = 1:le
    plot(selTr(i).points(:,1),selTr(i).points(:,2),'r-')
    plot(selTr(i).points(end,1),selTr(i).points(end,2),'r*')
end