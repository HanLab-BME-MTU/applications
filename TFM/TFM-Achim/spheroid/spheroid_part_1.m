function []=spheroid_part_1()
try
    load('fileAndFolderNames');
catch
    display('Browse to the project folder!')
    return;
end

load(path_displField);
load(path_forceField);

toDoList=1:length(forceField);

[maskSphrd,edgeSphrd]=spheroidDetect([],[],[],toDoList,[],[],1);

% cut off 3 grid spacings off the force field at each side of the field of
% view:

forceFieldSlim=forceField;
for iframe=toDoList
    minPos=min(forceField(iframe).pos,[],1);
    maxPos=max(forceField(iframe).pos,[],1);
    
    margin=3*forceField(iframe).par.gridSpacing;    
    cv=(forceField(iframe).pos(:,1)>minPos(1)+margin & forceField(iframe).pos(:,2)>minPos(2)+margin & forceField(iframe).pos(:,1)<maxPos(1)-margin & forceField(iframe).pos(:,2)<maxPos(2)-margin);
    
    forceFieldSlim(iframe).pos=forceField(iframe).pos(cv,:);
    forceFieldSlim(iframe).vec=forceField(iframe).vec(cv,:);
    
    forceFieldSlim(iframe).posShifted=forceField(iframe).posShifted(cv,:); %the original position counts here! Not the shifted one!
    forceFieldSlim(iframe).vecReIntp =forceField(iframe).vecReIntp(cv,:);
end

%**************************************************************************
% Plot the spheroid channel with segmentation and force field:
%**************************************************************************

inputFileList =getFileListFromFolder(['data',filesep,'6SpheroidFinal']);
target_dir_sphrdWithShiftedForces=[path_mechTFM,filesep,'sphrdWithShiftedForces'];

% Check if the dir exists:
if isdir(target_dir_sphrdWithShiftedForces)
    rmdir(target_dir_sphrdWithShiftedForces,'s');
    mkdir(target_dir_sphrdWithShiftedForces);
else
    mkdir(target_dir_sphrdWithShiftedForces);
end
plotCellsWithShiftedForces(inputFileList,forceFieldSlim,[],toDoList,2,edgeSphrd,target_dir_sphrdWithShiftedForces)


%**************************************************************************
% Plot the monolayer channel with segmentation and force field:
%**************************************************************************
inputFileList =getFileListFromFolder(['data',filesep,'6MonolayerFinal']);
target_dir_monolWithShiftedForces=[path_mechTFM,filesep,'monolWithShiftedForces'];

% Check if the dir exists:
if isdir(target_dir_monolWithShiftedForces)
    rmdir(target_dir_monolWithShiftedForces,'s');
    mkdir(target_dir_monolWithShiftedForces);
else
    mkdir(target_dir_monolWithShiftedForces);
end
plotCellsWithShiftedForces(inputFileList,forceFieldSlim,[],toDoList,2,edgeSphrd,target_dir_monolWithShiftedForces)

for iframe=toDoList;
    areaSphrd(iframe)= sum(maskSphrd(iframe).mat(:));
end
figure();
plot(areaSphrd,'ok')
xlabel('frame')
ylabel('spheroid spreading area [pix^2]')

display('all done!');