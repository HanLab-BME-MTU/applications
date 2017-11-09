% Given a set of arbitrary points in a movie and phyllopodia IDs, it'll draw segments
% that link the endpoints to the points

% The format of the variable originPoints should be like
% originPoints=[id1 X1 Y1; id2 X2 Y2; id3 X3 Y3]
% the IDs have to be real numbers in the fullData.dat file and the location
% obviously should be inside in the image

% Author: Antoine Godin
% godin.antoine@sympatico.ca

%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName=getFromINI(settingsRead, 'fileName');
frames=getFromINI(settingsRead, 'frames');
originPoints=getFromINI(settingsRead, 'idAndOrigin');

fullInfo=load([imagesFolder 'Tracks\' 'FullInfo.dat']);
segmentsList=[];
segments=[];

trackedFrames=unique(fullInfo(:,9)); % column 9 is the id

itwaitbar = waitbar(0,['Making Segments With Arbitrary Points']);     
for itFrame=1:frames
    thisFrame=imread([imagesFolder 'Skeletons\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    redCone=thisFrame(:,:,1);
    greenCone=thisFrame(:,:,2);
    blueCone=thisFrame(:,:,3);
    linesImage=zeros(size(thisFrame, 1), size(thisFrame, 2));
    if find(trackedFrames==itFrame)>0
        try
            thisFrameRows=fullInfo(find(fullInfo(:,9)==trackedFrames(itFrame)),:);
            col=0;
            for itArbitrary=1:size(originPoints, 1)
                if find(thisFrameRows(:,10)==originPoints(itArbitrary, 1))
                    col=col+1;
                    segments(col,:)=thisFrameRows(find(thisFrameRows(:,10)==originPoints(itArbitrary, 1)),:);
                    % replace the centroid columns for the arbitrary points
                    segments(col,6)=originPoints(itArbitrary, 2);
                    segments(col,7)=originPoints(itArbitrary, 3);
                end
            end
            for itSegment=1:size(segments, 1)
                linesImage=or(linesImage, func_Drawline(zeros(size(thisFrame, 1),...
                    size(thisFrame, 2)), segments(itSegment, 2),...
                segments(itSegment, 1), segments(itSegment, 7), ...
                segments(itSegment, 6), 1));
                linesImage=addNumbersImage(linesImage, segments(itSegment, 10), [7, 14],...
                    [segments(itSegment, 2), segments(itSegment, 1)]);
                segmentsList=[segmentsList; segments(itSegment, 10) ...
                    sqrt((segments(itSegment, 1)-segments(itSegment, 6))^2+ ...
                    (segments(itSegment, 2)-segments(itSegment, 7))^2) itFrame];
            end
            segments=[];
            redCone(linesImage~=0)=linesImage(linesImage~=0)*255;
            greenCone(linesImage~=0)=linesImage(linesImage~=0)*255;
        end
    end
    rbgCone(:,:,1)=redCone;
    rbgCone(:,:,2)=greenCone;
    rbgCone(:,:,3)=blueCone;
    imwrite(rbgCone, [imagesFolder '\Arbitrary\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'jpg')
    itwaitbar = waitbar(itFrame/frames,itwaitbar,['Making Segments With Arbitrary Points']);     
end
close(itwaitbar)

saveWithHeaders([imagesFolder 'Arbitrary\'], 'ArbitraySegments.dat', {'ID', 'Length',...
    'Frame'}, segmentsList);

%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'Arbitrary\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');