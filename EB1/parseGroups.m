function parseGroups(group,tracks)

if nargin == 0
    [fileName,dirName] = uigetfile('*.tif','Choose an image');
    load([dirName(1:end-8),'\point_files\config001_5p00_track_bidir.mat']);
    load([dirName(1:end-8),'\groups\group.mat']);
    traj=tracks(find([tracks.len]>=4));
end
leT = length(traj);
leG = length(group);

for i = 1:leG
    leGi = length([group(i).list]);
    for j = 1:(leGi-1)
        dy = traj(group(i).list(j+1)).points(1,1) - traj(group(i).list(j)).points(end,1);
        dx = traj(group(i).list(j+1)).points(1,2) - traj(group(i).list(j)).points(end,2);
        aVec = ([dx; dy])';
        magAvec = sqrt(sum(aVec.^2,2));

        dY = traj(group(i).list(j)).points(end,1) - traj(group(i).list(j)).points(1,1);
        dX = traj(group(i).list(j)).points(end,2) - traj(group(i).list(j)).points(1,2);
        traj_vec = [dX; dY];
        magTrajVec = sqrt(sum(traj_vec.^2,1));

        cos_ang = (aVec*traj_vec)./(magAvec*magTrajVec);
        if cos_ang > 0
            group(i).dire(j) = 1; % forward
        elseif cos_ang < 0
            group(i).dire(j) = -1; % backward
        end
    end
end




% calculate speeds
for i = 1 : le
    d = diff(tracks(i).points(:, 1:2));
    tracks(i).meanDisp = sum(sqrt(sum(d.^2, 2)))/(tracks(i).len-1) ;
end
