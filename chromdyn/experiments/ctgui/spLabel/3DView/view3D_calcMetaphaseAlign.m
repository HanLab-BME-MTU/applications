function view3D_calcMetaphaseAlign
%calculates rotation matrices to align mitotic spindles in view3D (at least 3 tags are needed!)

%load handles and data
labelPanelH = findall(0,'Tag','LabelPanel');
view3DH = GetUserData,'view3DGUIH');
view3D_handles = guidata(view3DH);
axisTags = view3D_handles.axisTagNumber;
cenTags = view3D_handles.nonAxisTagNumber;
axisCoordinates = view3D_handles.axisCoordinates;
rotMatrix = view3D_handles.rotMatrix;
goodTime = view3D_handles.goodTime;
nspots = view3D_handles.nspots;


idlist = GetUserData(labelPanelH,'idlist');

%if no more goodTime, no need to look further!
tmaxLoop = max(goodTime);
tmaxInit = size(idlist,2);

%init vars
cenCoord = zeros(tmaxInit,3);
e_perp = zeros(tmaxInit,3);
v_perp = zeros(tmaxInit,3);
e_3 = zeros(tmaxInit,3);
rotMatrixMeta = zeros(3,3,tmaxInit);

%calculate axis vectors
axisVectors = axisCoordinates(2,:,:)-axisCoordinates(1,:,:);
axisVectors = squeeze(axisVectors)'; %makes it a tmax by 3 matrix
[n_axisVectors,e_axisVectors] = normList(axisVectors);

%for each timepoint with at least three spots: calculate rotMatMeta, for all others, take rotMat from calcRotMat
for t = 1:tmaxLoop
    switch (nspots(t)>2)+(nspots(t)>0)
    case 2 %at least three spots
        idlist(t).linklist = sortrows(idlist(t).linklist,4);
        cenCoord(t,:) = mean(idlist(t).linklist(cenTags,9:11),1);
        %calculate vector through cen-center perpendicular to axis
        v_perp(t,:) = perpVector(idlist(t).linklist(axisTags(1),9:11),e_axisVectors(t,:),cenCoord(t,:));
        %norm perpendicular vector
        [n_perp(t,:),e_perp(t,:)] = normList(v_perp(t,:));
        %calculate third vector perpendicular to the two others
        v_3 = cross(e_axisVectors(t,:),e_perp(t,:));
        [dummy,e_3(t,:)] = normList(v_3);
        %rotation matrix is the inverse of the unit matrix specifying the coordinate system
        rotMatrixMeta(:,:,t) = [e_axisVectors(t,:);e_perp(t,:);e_3(t,:)];
    case 1 %1-2 spots: take from view3D_calcRotMat
        rotMatrixMeta(:,:,t) = rotMatrix(:,:,t);
    case 0 %no goodTime
    end
end

%store results
view3D_handles.rotMatrixMeta = rotMatrixMeta;
view3D_handles.cenCenterDistance = n_perp;

guidata(view3DH,view3D_handles);