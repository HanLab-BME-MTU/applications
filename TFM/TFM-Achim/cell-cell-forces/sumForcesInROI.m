function [sumForceVec netStruct]=sumForcesInROI(forceField,imageFileList,bwMask,toDoList,pixSize_mu)

%read in Stack of images:
if nargin < 2 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Ecad-Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
elseif isdir(imageFileList)
    imageFileList=getFileListFromFolder(imageFileList);
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin<4 || isempty(toDoList)
    toDoList=1:length(forceField);
end

if nargin<5 || isempty(pixSize_mu)
    pixSize_mu=0.163;
end

% Read image to retrieve crop coordinates
frame=1;
img=double(imread(imageFileList{1}));
h=figure;
set(h,'NumberTitle','off');
set(h,'Name','Please draw region of interest');
% Normalize img
img=(img-min(img(:)))/(max(img(:))-min(img(:)));

% Crop - if the user closes the window without drawing, roipoly will return an error
try
    [imgCropped,area]=imcrop(img);
catch
    uiwait(msgbox('No polygon selected. Quitting','Error','modal'));
    return
end
close(h);
% Round area (imcrop can give also non-integer boundaries)
area=round(area);
% Vertices
y0=area(2); y=area(2)+area(4);
x0=area(1); x=area(1)+area(3);

% Check boundaries
if y0<=0, y0=1; end
if x0<=0, x0=1; end
if y>=size(img,1), y=size(img,1); end
if x>=size(img,2), x=size(img,2); end

% Cut image
bwMask=zeros(size(img));
bwMask(y0:y,x0:x)=1;

frame=1;
figure()
imagesc(bwMask);
hold on
quiver(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1),forceField(frame).vec(:,2));
hold off;

area = sum(bwMask(:));
for frame=toDoList
    [currSumForceVec,~,method,~,~]=integrateForceField(forceField(frame).pos,forceField(frame).vec,bwMask,pixSize_mu);
    sumForceVec(frame).vec=currSumForceVec;
    sumForceVec(frame).vec_perpix=currSumForceVec/sum(bwMask(:));
    
    netStruct{frame}.stats.errorSumForce.vec        = currSumForceVec;
    netStruct{frame}.stats.errorSumForce.vec_perpix = currSumForceVec/area;
    netStruct{frame}.stats.errorSumForce.method     = method;
end
for frame=toDoList
    plot(sumForceVec(frame).vec(1),sumForceVec(frame).vec(2),'-*')
    hold on
end