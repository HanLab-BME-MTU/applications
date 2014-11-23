%% ---- parameters ------- %%

outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Supplemental Movies\';

useMov = true;
    
curvType = 'Max';

%% --- Curvature over time example ---- %%




movPath = 'W:\Hunter\orchestra_files_and_backup_merged\nih\20090605_LatestData_Blood\TdCAAX2\set3\pos4\movieData.mat';


MD = MovieData.load(movPath,0);

%%

if ~isempty(MD.eventTimes_)
    nFrames = MD.eventTimes_;
else
    nFrames = MD.nFrames_;
end

pixXY = MD.pixelSize_;
iChan = 1;

outFile = [outDir curvType ' curvature over time'];

frameRate = 3;
qual = 99;
if useMov
    MakeQTMovie('start',[outFile '.mov']);
    MakeQTMovie('quality',qual/100);
    MakeQTMovie('framerate',frameRate);
else
    writerObj = VideoWriter([outFile '.avi']);    
    writerObj.FrameRate= frameRate;
    writerObj.Quality = qual;
    open(writerObj)
end
    
currFig = fsFigure(.75);
cla;
currAx = get(currFig,'CurrentAxes');

imSize = MD.imSize_;
xl = [0 imSize(2)*pixXY/1e3];
yl = [0 imSize(1)*pixXY/1e3];
zl = [0 MD.nSlices_*pixXY/1e3];%mesh is made with isotropic voxels

camAng = 6.77118;
camTarg = [30.2285 29.52 5.14328];
camPos = [-195 138 356];





for j = 1:nFrames
    
   mg = MD.processes_{MD.getProcessIndex('MaskGeometry3DProcess',1,0)}.loadChannelOutput(iChan,j);
   mg.SmoothedSurface.vertices = mg.SmoothedSurface.vertices * pixXY / 1e3;   
   
      
   cla(currAx)
   
   switch curvType
       
       case 'Mean'
            cDat = mg.locAvgCurv.LocMeanMeanCurvature / pixXY * -1e3;            
            colormap(curvColormap)            
            cRange = [-.75 .75];
       case 'Max'
            cDat = mg.locAvgCurv.LocMeanMaxAbsCurvature / pixXY * 1e3;            
            colormap(jet(256))
            cRange = [.4 1.5];
   end
   
   
   pHan = patch(mg.SmoothedSurface,'EdgeColor','none','FaceColor','flat','FaceVertexCData',cDat);
   caxis(cRange)      
   
   set(pHan,'VertexNormals',mg.SurfaceNorms(:,[1 2 3]));
   axis equal
   set(currAx,'Projection','perspective')
   light
   l2 = light;
   lighting phong
   set(l2,'Position',[-1 0 1])
   
   xlim(xl),ylim(yl),zlim(zl)
   set(currFig,'color','w')
   set(currAx,'color','w')
   
   set(currAx,'CameraViewAngle',camAng)
   set(currAx,'CameraTarget',camTarg)
   set(currAx,'CameraPosition',camPos);
      
   if useMov
        MakeQTMovie('addfigure')
   else
        currFrame = getframe(currFig);
        writeVideo(writerObj,currFrame);
   end
   
    
end

if useMov
    MakeQTMovie('finish')
else
    close(writerObj)
end

%% ------- Curv category movie is in figures script because thats where the category code was ---- %%

