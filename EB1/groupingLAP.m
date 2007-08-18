function group = groupingLAP(load_mat)

% pass the TRACKS struct as input
nbFrames = 125;
LifeTime = 4; %4

MaxTimeSpan = 15;
coneAngPos = cos(pi/4);%45
coneAngNeg = cos(pi/18);%10
trackAng = cos(pi/3);%60
constVel = 1; %debug100
duConst = 1000000;
if nargin == 1
    [fileName,dirName] = uigetfile('*.mat','Choose a .mat file');
    load([dirName,filesep,fileName]);
else
    %     load(['/mnt/alex10/AlexData11/786O/786O_parental/786Opar_NaCl02_R3D/point_files/config001_5p00_track_bidir.mat']);
    load(['C:\amatov\data\786O\config001_5p00_track_bidir.mat']);
end

indx = find( [tracks.len] >= LifeTime);
TimeSpan = min(max([tracks.len]),MaxTimeSpan); 
traj = tracks(indx);
leIndx = length(indx);

% calculate trajectory information
for i = 1:leIndx
    traj(i).endID = traj(i).startID + traj(i).len - 1;
    dY = traj(i).points(end,1) - traj(i).points(1,1);
    dX = traj(i).points(end,2) - traj(i).points(1,2);
    traj(i).vec = [dX; dY];
    traj(i).vel = sqrt(dY^2+dX^2)/traj(i).len;
end
% exclude tracks with NO motion
traj = traj(find([traj.vel]>0)); % SOME DONT MOVE!?
leIndx = length(traj);

% get start & end points of all tracks for distance matrix calculation
for i = 1:leIndx
    E(i,:)=traj(i).points(end,1:2);
    B(i,:)=traj(i).points(1,1:2);
end

[listStartID,indxStartID] = sort([traj.startID]);
M = sparse(leIndx,leIndx); % angle multiplication matrix
C = sparse(leIndx,leIndx); % distance cut-off matrix

progressText(0,'First FOR-loop') % Create text
for i = 1:leIndx
    if TimeSpan > nbFrames-traj(i).endID
        T = nbFrames-traj(i).endID;
    else
        T = TimeSpan;
    end
    if (traj(i).endID < (nbFrames - LifeTime))

        aux1 = find(listStartID>(traj(i).endID+1));
        firstIndx =  aux1(1);

        aux2 = find(listStartID>(traj(i).endID+T+1));
        if ~isempty(aux2)
            finalIndx =  aux2(1)-1;
        else
            finalIndx = length(listStartID);
        end

        a = indxStartID(firstIndx:finalIndx);

        magTrajVec = sqrt(sum(traj(i).vec.^2,1));
        magTrajA = sqrt(sum([traj(a).vec]'.^2,2));
        cos_trackAng = ([traj(a).vec]'*traj(i).vec)./(magTrajA*magTrajVec);
        
        listA = find(cos_trackAng>trackAng);
        aa = a(listA);
        cos_trA = cos_trackAng(listA);
   
        if ~isempty(aa)
            dxCo = []; dyCo = [];
            for j = 1:length(aa)
                dyCo(j) = traj(aa(j)).points(1,1) - traj(i).points(end,1); % get it out of the loop - to the other loop up
                dxCo(j) = traj(aa(j)).points(1,2) - traj(i).points(end,2);
            end
            aVecs = ([dxCo; dyCo])';

            magAvecs = sqrt(sum(aVecs.^2,2));
            magAvecs(find(magAvecs==0)) = 0.001; % some new/t+4 traj begin where an old/t one ended
            cos_coneAng = (aVecs*traj(i).vec)./(magAvecs*magTrajVec);

            posC = find(cos_coneAng>coneAngPos);
            negC = find(cos_coneAng<-coneAngNeg);

            if ~isempty(posC)
                C(i,aa(posC)) = traj(i).vel * sqrt([traj(aa(posC)).startID]-traj(i).endID) ;%* constVel;
                M(i,aa(posC)) = abs(cos_trA(posC) - cos_coneAng(posC));
            end
            if ~isempty(negC)
                C(i,aa(negC)) = traj(i).vel * sqrt([traj(aa(negC)).startID]-traj(i).endID) ;%* constVel;
                M(i,aa(negC)) = abs(cos_trA(negC) - abs(cos_coneAng(negC)));
            end
        end
    end
    progressText(i/leIndx);
end

R = max(C(:));
D=createSparseDistanceMatrix(E,B,R); % dimentions - leIndx
radIndx = find(D<=C);
C(radIndx) = D(radIndx).*M(radIndx); % overwrite old cut-off matrix to generate cost matrix

[links12, links21] = lap(C,-1,0,1);
lnk = links12(1:leIndx);

lnkIndx = find(lnk<=leIndx);
leLnkIndx = length(lnkIndx);
group=struct('list',[]);
assocTracks = zeros(leIndx,1);
grNb = 0; 

progressText(0,'Second FOR-loop') % Create text
for i = 1:leLnkIndx
    k = lnkIndx(i);
    if assocTracks(k) == 0
        if assocTracks(lnk(k)) == 0
            grNb = grNb + 1;
            assocTracks(k) = grNb;
            assocTracks(lnk(k)) = grNb;
            group(grNb).list = [k lnk(k)];
            traj(k).next = lnk(k);
            traj(lnk(k)).prev = k;
            traj(k).g_nb = grNb;
            traj(lnk(k)).g_nb = grNb;
        else
            assocTracks(k) = assocTracks(lnk(k));
            group(assocTracks(lnk(k))).list = [k,group(assocTracks(lnk(k))).list];
            traj(k).next = lnk(k);
            traj(lnk(k)).prev = k;
            traj(k).g_nb = assocTracks(lnk(k));
        end
    else
        assocTracks(lnk(k)) = assocTracks(k);
        group(assocTracks(k)).list = [group(assocTracks(k)).list,lnk(k)];
        traj(k).next = lnk(k);
        traj(lnk(k)).prev = k;
        traj(lnk(k)).g_nb = assocTracks(k);
    end
    progressText(i/leLnkIndx);
end

grNb = length(group);
figure
for i = 1:grNb
    k = group(i).list;
    le_gr = length(k);
    plot(traj(k(1)).points(1,2),traj(k(1)).points(1,1),'ks') 
    for j = 1:le_gr
        plot(traj(k(j)).points(:,2),traj(k(j)).points(:,1),'r-')
        hold on
        plot(traj(k(j)).points(end,2),traj(k(j)).points(end,1),'r*')
        plot([traj(k(j)).points(1,2),traj(k(j)).points(end,2)],[traj(k(j)).points(1,1),traj(k(j)).points(end,1)],'k:')
        if j < le_gr
            plot([traj(k(j)).points(end,2),traj(k(j+1)).points(1,2)],[traj(k(j)).points(end,1),traj(k(j+1)).points(1,1)],'g-')
            plot(traj(k(j+1)).points(1,2),traj(k(j+1)).points(1,1),'g*')
        end
    end
end
hold off
disp(sprintf('Number of groups %d',grNb));