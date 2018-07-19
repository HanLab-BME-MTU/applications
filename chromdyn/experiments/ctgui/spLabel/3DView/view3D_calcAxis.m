function view3D_calcAxis(updateSwitch)
%using idlist from labelpanel and axis definition from view3D, view3D_calcAxis saves start and end coordinates for axis

if nargin==0|isempty(updateSwitch)
    updateSwitch = 0;
end

%load handles/data
labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
view3DH = GetUserData(labelPanelH,'view3DGUIH');
view3D_handles = guidata(view3DH);
axisTags = view3D_handles.axisTagNumber;
idlist = GetUserData(labelPanelH,'idlist');

tmax = size(idlist,2);

%centering
centerMode = get(view3D_handles.view3D_center_PD,'Value');
switch centerMode
    case 1 %center none
        centerString = '';
    case 2 %center centroid
        centerString = '-idlist(t).centroid';
    case 3 %center tag
        tagNum = view3D_handles.centerTagNumber; % = log2(tagColor)+1
        centerString = ['-idlist(t).linklist(',num2str(tagNum),',9:11)'];
end

if updateSwitch %if handles are being updated, update line coordinates according to align status, too
    alignMode = get(view3D_handles.view3D_align_PD,'Value');
    switch alignMode
        case 1
            alignString = '';
            transposeString = '';
        case 2 %align axis
            rotMatrix = view3D_handles.rotMatrix;
            alignString = 'rotMatrix(:,:,t)*';
            transposeString = '''';
        case 3 %align metaphase
            rotMatrixMeta = view3D_handles.rotMatrixMeta;
            alignString = 'rotMatrixMeta(:,:,t)*';
            transposeString = ''''; 
    end
else
    alignString = '';
    transposeString = '';
end

%build axis-assignment string
coordString1 = ['axisCoordinates(1,:,t) = ',alignString,'(idlist(t).linklist(axisTags(1),9:11)',centerString,')',transposeString,'; '];
coordString2 = ['axisCoordinates(2,:,t) = ',alignString,'(idlist(t).linklist(axisTags(2),9:11)',centerString,')',transposeString,'; '];


%sort linklists according to tag color and read axis coordinates 
axisCoordinates = zeros(2,3,tmax);

for t = 1:tmax
    if ~isempty(idlist(t).linklist)
        idlist(t).linklist = sortrows(idlist(t).linklist,4);
        eval(coordString1);
        eval(coordString2)
    end
end

view3D_handles.axisCoordinates = axisCoordinates;
view3D_handles.axisDefined = 1;

guidata(view3DH,view3D_handles);