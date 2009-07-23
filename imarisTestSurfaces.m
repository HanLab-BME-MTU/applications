function imarisTestSurfaces(obj)


if nargin < 1
    obj = [];
end


%% Create surfaces and show in Matlab

% make a coordinate system
[xx,yy,zz] = meshgrid(1:100,1:100,1:100);

% create a volume (xyztc)
volume = zeros(100,100,100,1,2);
% fill the first frame with a Gaussian
volume(:,:,:,1,1) = exp( -(xx-50).^2/100 -(yy-50).^2/400 -(zz-50).^2/200);
% fill the second frame with a cuboid
volume(30:50,23:78,40:60,1,2) = 1;

% create surfaces, plot
for t=1:2
    [isoData(t).F isoData(t).V] = isosurface(xx,yy,zz,volume(:,:,:,1,t),0.33);
    isoData(t).nV = size(isoData(t).V,1);
    isoData(t).nF = size(isoData(t).F,1);
    isoData(t).T = t;
    isoData(t).N = isonormals(volume(:,:,:,1,t),isoData(t).V);
    
    figure('Name',sprintf('frame %i',t))
    ph = patch('Faces',isoData(t).F,'Vertices',isoData(t).V);
    isonormals(volume(:,:,:,1,t),isoData(t).V);
    set(ph,'EdgeColor','none','FaceColor','r','FaceLighting','phong')
    camlight, lighting gouraud
    xlim([1,100])
    ylim([1,100])
    zlim([1,100])
end

% catenate data to be used for addSurfaceList
imaData.V = cat(1,isoData(:).V);
imaData.F = cat(1,isoData(:).F);
imaData.N = cat(1,isoData(:).N);
imaData.nF = cat(1,isoData(:).nF);
imaData.nV = cat(1,isoData(:).nV);
imaData.T = cat(1,isoData(:).T);

dataSize = ones(1,5);
dataSize(1:ndims(volume)) = size(volume);

%% Start Imaris and show volume
if isempty(obj)
    imaApp = imarisStartNew;
    % % create dataSet in Imaris
    imaDataSet = imaApp.mFactory.CreateDataSet;
    % careful: imaris gets data as xyzct
    imaDataSet.Create('eTypeFloat',dataSize(1),dataSize(2),dataSize(3),dataSize(4),dataSize(5));
    imaDataSet.SetData(single(volume));
    % make a surpass scene
    imaSurpassScene = imaApp.mFactory.CreateDataContainer;
    % add dataSet to Imaris
    imaApp.mDataSet = imaDataSet;
else
    imaApp = obj.imageData.showImaris;
end

% fill surpass scene with light and frame
imaLight = imaApp.mFactory.CreateLightSource;
imaSurpassScene.AddChild(imaLight);
imaFrame = imaApp.mFactory.CreateFrame;
imaSurpassScene.AddChild(imaFrame);
% add volume
imaVolume = imaApp.mFactory.CreateVolume;
imaSurpassScene.AddChild(imaVolume);

% add surpass scene and set view
imaApp.mSurpassScene = imaSurpassScene;
imaApp.mViewer = 'eViewerSurpass';

%% Show surface in Imaris
% get surpass scene from Imaris
imaSurpassScene = imaApp.mSurpassScene;
% add surfaces
imaSurfaces = imaApp.mFactory.CreateSurfaces;
for t=1:2
    imaSurfaces.AddSurface(isoData(t).V(:,[2,1,3]),isoData(t).F-1,isoData(t).N(:,[2,1,3]),isoData(t).T-1);
end
imaSurpassScene.AddChild(imaSurfaces);


%% subfunction ImarisStartNew
function imarisApplication = imarisStartNew(assigninBase)
%IMARISSTARTNEW starts a new Imaris and stores the handle in the workspace
% naming it imarisApplication, or imarisApplication#, where # is a number
%
% SYNOPSIS imarisApplication = imarisStartNew
%
% INPUT    assigninBase (opt): true if Imaris should assign handle in base
%                              Default: true
%
% OUTPUT   imarisApplication : Handle to the new imaris session
%
%c: jonas 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 || isempty(assigninBase)
    assigninBase = true;
end
if ~ispc
    error('Imaris only runs under Windows!')
end

imarisApplication = actxserver('Imaris.Application');
imaAppName = 'imarisApplication';
% make sure we do not accidentially overwrite an imaris application
if evalin('base','exist(''imarisApplication'',''var'')')
    num = 2;
    while evalin('base',['exist(''imarisApplication',num2str(num),''',''var'')'])
        num = num+1;
    end
    imaAppName = [imaAppName, num2str(num)];
end
imarisApplication.mVisible = 1;
if assigninBase
    assignin('base',imaAppName,imarisApplication);
    disp(sprintf(['The handle to the current imaris is ''%s''\n',...
        'Delete it to close this session of Imaris'],imaAppName));
end

