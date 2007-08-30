function selectTracks(axis,traj,I,n)

if nargin == 0
    [fileName,dirName] = uigetfile('*.tif','Choose an image');
    I = imread([dirName,fileName]);
    load([dirName(1:end-8),'\point_files\config001_6p00_track_bidir.mat']);
    load([dirName(1:end-8),'\poles\axis.mat']);
    traj=tracks(find([tracks.len]>=2));
    n = 20; LTw = 7;
    leTraj = length(traj);
    for i = 1:leTraj
        traj(i).endID = traj(i).startID + traj(i).len - 1;
        dY = traj(i).points(end,1) - traj(i).points(1,1);
        dX = traj(i).points(end,2) - traj(i).points(1,2);
        traj(i).vec = [dX; dY];
        traj(i).vel = sqrt(dY^2+dX^2)/traj(i).len;
    end
    % exclude tracks with NO motion
    traj = traj(find([traj.vel]>0)); % SOME DONT MOVE!?
    % get the tracks only for the window
    indxT = find([traj.endID]>=(n+LTw) & [traj.startID]<(n+LTw+3));
    traj = traj(indxT);
end
leTraj = length(traj);
% calculate poles information for the n-th set of pole coordinates
dYp = poleAxis(n,3) - poleAxis(n,1);
dXp = poleAxis(n,4) - poleAxis(n,2);
poleVec = [dXp; dYp];
for i = 1:leTraj
    magTrajVec = sqrt(sum(traj(i).vec.^2,1));
    % cluster tracks based on distance to each pole of track's beginning
    dy1 = traj(i).points(1,1) - poleAxis(n,1);
    dx1 = traj(i).points(1,2) - poleAxis(n,2);
    aVecs1 = ([dx1; dy1])';
    magAvecs1 = sqrt(sum(aVecs1.^2,2));
    cos_clustAng1 = (aVecs1*traj(i).vec)./(magAvecs1*magTrajVec);

    dy2 = traj(i).points(1,1) - poleAxis(n,3);
    dx2 = traj(i).points(1,2) - poleAxis(n,4);
    aVecs2 = ([dx2; dy2])';
    magAvecs2 = sqrt(sum(aVecs2.^2,2));
    cos_clustAng2 = (aVecs2*traj(i).vec)./(magAvecs2*magTrajVec);

    if cos_clustAng1 < 0 & cos_clustAng2 < 0
        traj(i).pol = 0;
    elseif cos_clustAng1 > 0 & (cos_clustAng2) < 0
        traj(i).pol = 1;
    elseif cos_clustAng2 > 0 & (cos_clustAng1) < 0
        traj(i).pol = 2;
    elseif sqrt(dy1^2+dx1^2) > (sqrt(dy2^2+dx2^2))
        traj(i).pol = 2;
    elseif sqrt(dy2^2+dx2^2) > (sqrt(dy1^2+dx1^2))
        traj(i).pol = 1;
    else
        traj(i).pol = 3;
    end
end
indx1 = find([traj.pol]==1);
nbTrPole1 = length(indx1);
indx2 = find([traj.pol]==2);
nbTrPole2 = length(indx2);
indx0 = find([traj.pol]==3);
nbTrPole0 = length(indx0)
ALL_TRACKS = nbTrPole1 + nbTrPole2
% plot the poles and the tracks
figure,imshow(I,[])
hold on
plot(poleAxis(n,1),poleAxis(n,2),'gs')
plot(poleAxis(n,3),poleAxis(n,4),'gs')
for i = 1:nbTrPole0
    plot(traj(indx0(i)).points(:,1),traj(indx0(i)).points(:,2),'y-')
    plot(traj(indx0(i)).points(end,1),traj(indx0(i)).points(end,2),'y*')
end
for i = 1:nbTrPole1
    plot(traj(indx1(i)).points(:,1),traj(indx1(i)).points(:,2),'b-')
    plot(traj(indx1(i)).points(end,1),traj(indx1(i)).points(end,2),'b*')
end
for i = 1:nbTrPole2
    plot(traj(indx2(i)).points(:,1),traj(indx2(i)).points(:,2),'r-')
    plot(traj(indx2(i)).points(end,1),traj(indx2(i)).points(end,2),'r*')
end
figure,imshow(I,[])
hold on
TRACKS_OUT = 0;
%-------------------------------------------------------------------------
traj1 = traj(find([traj.pol]==1)); % look at pole1
% calculate the dot product
magPoleVec = sqrt(sum(poleVec.^2,1));
magTrajVec = sqrt(sum([traj1.vec]'.^2,2));
cos_trackAng1 = (poleVec'*[traj1.vec])./(magPoleVec*magTrajVec)';
list1 = find(cos_trackAng1<0); % for pole1
list2 = find(cos_trackAng1>=0);
list_out = list1;
list_in = list2;
% plot the poles and the tracks
plot(poleAxis(n,1),poleAxis(n,2),'gs')
plot(poleAxis(n,3),poleAxis(n,4),'gs')
for i = 1:length(list1)
    plot(traj1(list1(i)).points(:,1),traj1(list1(i)).points(:,2),'r-')
    plot(traj1(list1(i)).points(end,1),traj1(list1(i)).points(end,2),'r*')
end
%--------------------------------------------------------------------------
traj2 = traj(find([traj.pol]==2)); % look at pole2
% calculate the dot product
magPoleVec = sqrt(sum(poleVec.^2,1));
magTrajVec = sqrt(sum([traj2.vec]'.^2,2));
cos_trackAng2 = (poleVec'*[traj2.vec])./(magPoleVec*magTrajVec)';
list12 = find(cos_trackAng2>0); % for pole2
list22 = find(cos_trackAng2<=0);
list_out = [list_out,list12];
list_in = [list_in,list22];
for i = 1:length(list12)
    plot(traj2(list12(i)).points(:,1),traj2(list12(i)).points(:,2),'r-')
    plot(traj2(list12(i)).points(end,1),traj2(list12(i)).points(end,2),'r*')
end
%--------------------------------------------------------------------------
figure,imshow(I,[])
hold on
for i = 1:length(list2)
    plot(traj1(list2(i)).points(:,1),traj1(list2(i)).points(:,2),'r-')
    plot(traj1(list2(i)).points(end,1),traj1(list2(i)).points(end,2),'r*')
end
for i = 1:length(list22)
    plot(traj2(list22(i)).points(:,1),traj2(list22(i)).points(:,2),'r-')
    plot(traj2(list22(i)).points(end,1),traj2(list22(i)).points(end,2),'r*')
end
%--------------------------------------------------------------------------
NB_TRACKS_OUT = length(list1) + length(list12)
NB_TRACKS_IN = length(list2) + length(list22)
TRACKS_IN_PERC = NB_TRACKS_IN/(NB_TRACKS_IN+NB_TRACKS_OUT)
traj_out = [traj1(list1) traj2(list12)];
MEAN_LENGTH_OUT = mean([traj_out.len])
traj_in = [traj1(list2) traj2(list22)];
MEAN_LENGTH_IN = mean([traj_in.len])
RATION_LE_IN_TO_OUT = MEAN_LENGTH_IN/MEAN_LENGTH_OUT

traj



% selTr = traj(list1);
% le = length(list1);
% TRACKS_OUT = le + TRACKS_OUT;
% 
% % plot the poles and the tracks
% plot(poleAxis(n,1),poleAxis(n,2),'gs')
% plot(poleAxis(n,3),poleAxis(n,4),'gs')
% for i = 1:le
%     plot(selTr(i).points(:,1),selTr(i).points(:,2),'r-')
%     plot(selTr(i).points(end,1),selTr(i).points(end,2),'r*')
% end
% 
% le_in = length(list_in)
% figure,imshow(I,[])
% hold on
% plot(poleAxis(n,1),poleAxis(n,2),'gs')
% plot(poleAxis(n,3),poleAxis(n,4),'gs')
% for i = 1:le_in
%     plot(traj(list_in(i)).points(:,1),traj(list_in(i)).points(:,2),'b-')
%     plot(traj(list_in(i)).points(end,1),traj(list_in(i)).points(end,2),'b*')
% end
% TRACKS_OUT
% PERC_TRACKS_OUT = TRACKS_OUT / ALL_TRACKS
% MEAN_LENGTH = mean([traj.len])
% MEAN_LENGHT_OUT = mean([traj(list_out).len])
% MEAN_LENGHT_IN = mean([traj(list_in).len])