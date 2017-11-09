%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName     = getFromINI(settingsRead, 'fileName');
frames       = getFromINI(settingsRead, 'frames');
originPoints = getFromINI(settingsRead, 'idAndOrigin');
        % put origin and joints in the propper form
        xNear = originPoints(:,2);
        yNear = originPoints(:,3);        

%% Load data
skeletonEnds      = load([imagesFolder 'Tracks\' fileName 'Ends.dat']);
skeletonJoints    = load([imagesFolder 'Tracks\' fileName 'Joints.dat']);
fullInfo          = load([imagesFolder 'Tracks\' 'FullInfo.dat']);
trackedFrames     = unique(fullInfo(:,9));
segmentsList      = [];
phyllopodiaLength = [];

trackedFrames     = unique(fullInfo(:,9)); % column 9 is the id

warning off
itwaitbar = waitbar(0,['Making Segments by following Skeleton']);     
for itFrame=1:frames
    [itFrame,frames]
    endsToTrack           = [];
    thisPhyllopodiaLength = [];   
    thisSkeletons = imread([imagesFolder 'nakedSkeletons\' fileName num2str(itFrame-1, '%04.f') '.tif']);
    rgbPath         = zeros(size(thisSkeletons,1),size(thisSkeletons,2),3);
    rgbPath(:,:,1)  = double(imread([imagesFolder fileName num2str(itFrame-1, '%04.f') '.tif']))/255;
    try
        thisFrameRows = fullInfo(find(fullInfo(:,9)==trackedFrames(itFrame)),:);
        % find origin coordinates
        counter = 1;
        for itOrigin = 1:size(originPoints, 1)
            if size(find(thisFrameRows(:,10) == originPoints(itOrigin, 1)),1)>0
                fila = thisFrameRows(find(thisFrameRows(:,10) == originPoints(itOrigin, 1)),:);
                endsToTrack(counter, 1) = fila(1,1);
                endsToTrack(counter, 2) = fila(1,2);
                endsToTrack(counter, 3) = fila(1,10);
                counter = counter + 1;
            end
        end    
        if size(endsToTrack,1)>1 %only does the tracing if there are ends to track in this frame
            thisFrameEnds   = [skeletonEnds(find(skeletonEnds(:,3)==itFrame),1),skeletonEnds(find(skeletonEnds(:,3)==itFrame),2)];
            thisSkeletonsWithoutJoints = thisSkeletons;
            [connectivityMatrix,russianForm,PositionVectices] = GetConnectivityPointbyPoint(thisSkeletons);
            closestPoints = getClosestJoint(xNear,yNear,PositionVectices);      
            for itPhyllopodia=1:size(endsToTrack, 1)
                idBeg = closestPoints(find(endsToTrack(itPhyllopodia,3) == originPoints(:,1)));
                idEnd = find(and(PositionVectices(:,1)==endsToTrack(itPhyllopodia,1),PositionVectices(:,2)==endsToTrack(itPhyllopodia,2)));                  
                X1t = max([min([PositionVectices(idBeg(1),2),PositionVectices(idEnd(1),2)])-5,1]);
                X2t = min([max([PositionVectices(idBeg(1),2),PositionVectices(idEnd(1),2)])+5,size(thisSkeletons,1)]);
                Y1t = max([min([PositionVectices(idBeg(1),1),PositionVectices(idEnd(1),1)])-5,1]);
                Y2t = min([max([PositionVectices(idBeg(1),1),PositionVectices(idEnd(1),1)])+5,size(thisSkeletons,2)]);          
                    Connected = 0;  
                    while Connected == 0
                        thisSkeletons_T = thisSkeletons([X1t:X2t],[Y1t:Y2t]);
                        try
                            [connectivityMatrix_T,russianForm_T,PositionVectices_T] = GetConnectivityPointbyPoint(thisSkeletons_T);                    
                            idBeg_T = find(and(PositionVectices_T(:,1)==PositionVectices(idBeg(1),1)-Y1t+1,PositionVectices_T(:,2)==PositionVectices(idBeg(1),2)-X1t+1));                          
                            idEnd_T = find(and(PositionVectices_T(:,1)==endsToTrack(itPhyllopodia,1)-Y1t+1,PositionVectices_T(:,2)==endsToTrack(itPhyllopodia,2)-X1t+1));       
                            [dSP_T,shotestPath_T] = grShortPath(russianForm_T,idBeg_T,idEnd_T);                
                            if length(shotestPath_T) == 0
                                X1t = max([1,round(X1t)-5]);
                                X2t = min([round(X2t)+5,size(thisSkeletons,1)]);
                                Y1t = max([1,round(Y1t)-5]);
                                Y2t = min([round(Y2t)+5,size(thisSkeletons,2)]);
                            else
                                Connected = 1;
                            end
                        catch
                            X1t = max([1,round(X1t)-5]);
                            X2t = min([round(X2t)+5,size(thisSkeletons,1)]);
                            Y1t = max([1,round(Y1t)-5]);
                            Y2t = min([round(Y2t)+5,size(thisSkeletons,2)]);
                        end
                    end
                shortWayImage = zeros(size(thisSkeletons));
                for itSegments=1:size(shotestPath_T,2)
                    shortWayImage(PositionVectices_T(shotestPath_T(itSegments),2)+X1t,PositionVectices_T(shotestPath_T(itSegments),1)+Y1t) = 1;
                end            
                thisPhyllopodiaLength(itPhyllopodia,1) = getDistanceFollowingSkeleton(PositionVectices_T,shotestPath_T);
                thisPhyllopodiaLength(itPhyllopodia,2) = endsToTrack(itPhyllopodia,3);
                thisPhyllopodiaLength(itPhyllopodia,3) = itFrame;
                rgbPath(:,:,2) =  or(rgbPath(:,:,2),shortWayImage);
            end        
            phyllopodiaLength = [phyllopodiaLength; thisPhyllopodiaLength];
        end
        rgbPath(:,:,3) = thisSkeletons;
    end
    imwrite(rgbPath, [imagesFolder 'FollowSkeleton\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'JPG')
    itwaitbar = waitbar(itFrame/frames,itwaitbar,['Making Segments by following Skeleton']);
end
close(itwaitbar)
saveWithHeaders([imagesFolder 'FollowSkeleton\'], 'LengthAlongPhyllopodia.dat', {'Length', 'Id','Frame'}, phyllopodiaLength);

%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'FollowSkeleton\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');