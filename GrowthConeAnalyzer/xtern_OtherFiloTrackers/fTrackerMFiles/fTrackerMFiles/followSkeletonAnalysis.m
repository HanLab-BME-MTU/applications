%% Read parameters form settings file

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName=getFromINI(settingsRead, 'fileName');
frames=getFromINI(settingsRead, 'frames');
originPoints=getFromINI(settingsRead, 'idAndOrigin');

%% Load data
skeletonEnds=load([imagesFolder 'Tracks\' fileName 'Ends.dat']);
skeletonJoints=load([imagesFolder 'Tracks\' fileName 'Joints.dat']);
fullInfo=load([imagesFolder 'Tracks\' 'FullInfo.dat']);
trackedFrames=unique(fullInfo(:,9));
segmentsList=[];
phyllopodiaLength=[];

trackedFrames=unique(fullInfo(:,9)); % column 9 is the id

warning off
itwaitbar = waitbar(0,['Making Segments by following Skeleton']);     
for itFrame=1:frames
    endsToTrack=[];
    thisPhyllopodiaLength=[];
    
    thisSkeletons=imread([imagesFolder 'nakedSkeletons\' fileName num2str(itFrame-1, '%04.f') '.tif']);
    thisFrameRows=fullInfo(find(fullInfo(:,9)==trackedFrames(itFrame)),:);

    % find origin coordinates

    counter=1;
    for itOrigin=1:size(originPoints, 1)
        if size(find(thisFrameRows(:,10)==originPoints(itOrigin, 1)),1)>0
            fila=thisFrameRows(find(thisFrameRows(:,10)==originPoints(itOrigin, 1)),:);
            endsToTrack(counter, 1)=fila(1,1);
            endsToTrack(counter, 2)=fila(1,2);
            endsToTrack(counter, 3)=fila(1,10);
            counter=counter+1;
        end
    end
    rgbPath=zeros(size(thisSkeletons,1),size(thisSkeletons,2),3);
    rgbPath(:,:,1)=double(imread([imagesFolder fileName num2str(itFrame-1, '%04.f') '.tif']))/255;
    
    
    if size(endsToTrack,1)>0 %only does the tracing if there are ends to track in this frame
        % Ends and Joints
        thisFrameEnds=[skeletonEnds(find(skeletonEnds(:,3)==itFrame),1),...
            skeletonEnds(find(skeletonEnds(:,3)==itFrame),2)];
        thisFrameJoints=[skeletonJoints(find(skeletonJoints(:,3)==itFrame),1),...
            skeletonJoints(find(skeletonJoints(:,3)==itFrame),2)];
        thisSkeletonsWithoutJoints=thisSkeletons;

        % put origin and joints in the propper form
        for itOrigins=1:size(originPoints, 1)
            xNear(itOrigins)=originPoints(itOrigins, 2);
            yNear(itOrigins)=originPoints(itOrigins, 3);
        end

        closestJoints=getClosestJoint(xNear, yNear, thisFrameJoints);

        vertexesAndNeighbours={};

        for itJoints=1:size(thisFrameJoints, 1)
            thisSkeletonsWithoutJoints(thisFrameJoints(itJoints,2)-1:thisFrameJoints(itJoints,2)+1,...
                thisFrameJoints(itJoints,1)-1:thisFrameJoints(itJoints,1)+1)=0;
        end

        for itJoints=1:size(thisFrameJoints, 1)
            vertexesAndNeighbours{size(vertexesAndNeighbours,2)+1}=...
                [thisFrameJoints(itJoints, 1) thisFrameJoints(itJoints, 2);...
                getNeighbours([thisFrameJoints(itJoints, 1)...
                thisFrameJoints(itJoints, 2)], thisSkeletonsWithoutJoints)];
        end

        for itEnds=1:size(thisFrameEnds, 1)
                vertexesAndNeighbours{size(vertexesAndNeighbours,2)+1}=...
                    [thisFrameEnds(itEnds, 1) thisFrameEnds(itEnds, 2); ...
                    thisFrameEnds(itEnds, 1) thisFrameEnds(itEnds, 2)];
        end

        [labeledSkeleton, numbers]=bwlabel(thisSkeletonsWithoutJoints, 8);
    %    figure(1);imshow(label2rgb(labeledSkeleton, 'hsv', 'k', 'shuffle'));

        [connectivityMatrix, idConeexions]=getConnectivityAndColor(vertexesAndNeighbours, labeledSkeleton);
%         [connectivityMatrix,russianForm,PositionVectices] = GetConnectivityPointbyPoint(matrix)

        russianForm=getRussianForm(sparse(connectivityMatrix));

        %Add the two pixels I removed
        
        % When Chris reported the bugs I found this
        % russianForm(:,3)=russianForm(:,2)+2;. I will change it

        russianForm(:,3)=russianForm(:,3)+2;
       
        rgbPath(:,:,3)=thisSkeletons;
        
        for itPhyllopodia=1:size(endsToTrack, 1)
            idEnd=find(and(thisFrameEnds(:,1)==endsToTrack(itPhyllopodia,1),...
                thisFrameEnds(:,2)==endsToTrack(itPhyllopodia,2)));
            
            [dSP,shotestPath]=grShortPath(russianForm,idEnd+...
                size(thisFrameJoints, 1),closestJoints(itPhyllopodia));
            
            shortWayImage=zeros(size(labeledSkeleton));

            for itSegments=1:size(shotestPath, 2)-1
                shortWayImage(labeledSkeleton==idConeexions(shotestPath(itSegments), shotestPath(itSegments+1)))=1;
            end

            thisPhyllopodiaLength(itPhyllopodia,1)=sum(sum(shortWayImage));
            thisPhyllopodiaLength(itPhyllopodia,2)=endsToTrack(itPhyllopodia,3);
            thisPhyllopodiaLength(itPhyllopodia,3)=itFrame;
                       
            rgbPath(:,:,2)=or(rgbPath(:,:,2), shortWayImage);
        end
        
        phyllopodiaLength=[phyllopodiaLength; thisPhyllopodiaLength];
    
        rgbPath(:,:,2)=rgbPath(:,:,2);
        imwrite(rgbPath, [imagesFolder 'FollowSkeleton\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'JPG')
    else % no Ends to Track
        rgbPath=zeros(size(labeledSkeleton,1),size(labeledSkeleton,2),3);
        rgbPath(:,:,2)=thisSkeletons;
        imwrite(rgbPath, [imagesFolder 'FollowSkeleton\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'JPG')
    end
    close(itwaitbar)
    itwaitbar = waitbar(itFrame/frames,['Making Segments by following Skeleton']);
end
close(itwaitbar)
saveWithHeaders([imagesFolder 'FollowSkeleton\'], 'LengthAlongPhyllopodia.dat', {'Length', 'Id',...
    'Frame'}, phyllopodiaLength);

%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'FollowSkeleton\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');
