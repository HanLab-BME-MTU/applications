function [stateIdx,groupIdxArray] = trajectoryAnalysisMainCalcStatsFindGroups(dataListG,state)
%TAMCS-FINDGROUPS finds growth and shrinkage groups in dataLists. 
%Pauses of any length and undetermined periods with a length of one
%interval within a group are included

stateIdx = find(dataListG(:,3)==state);

%state belonging to the same phase are either direcly one after the other
%or separated by a pause -> make list of ones and zeros and find the groups


%try to find groups
gGroups = findGroups(dataListG(:,3),1);

%search where the groups are separated by just one other unit
%then check whether these are pauses or single undetermined
%then look through the list and fuse entries

gGroupsDiff  = gGroups(2:end,2)-gGroups(1:end-1,1);
gGroupsDiff2 = find(gGroupsDiff==2);

%lookfor pauses (gpgStates points to gGroupsDiff2, which points to gGroups, which points to dataListG)
gpgStates = find(dataListG(gGroups(gGroupsDiff2,2)+1,3)==3 | dataListG(gGroups(gGroupsDiff2,2)+1,3)==0);

%take again list for stateIdx
stateListState = dataListG(:,3);
%change pauses to state
stateListState(gGroups(gGroupsDiff2(gpgStates),2)+1) = state;

%find groups again
groupIdxArray = findGroups(stateListState,state);
