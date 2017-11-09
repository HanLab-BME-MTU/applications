% Saves a movie with segment that go from the center of mass to the
% skeleton ends

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName=getFromINI(settingsRead, 'fileName');
frames=getFromINI(settingsRead, 'frames');

fullInfo=load([imagesFolder 'Tracks\' 'FullInfo.dat']);

trackedFrames=unique(fullInfo(:,9));


%% Make segment images
itwaitbar = waitbar(1/2,'Making segment images');

for itFrame=1:frames
    thisFrame=imread([imagesFolder 'Skeletons\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    redCone=thisFrame(:,:,1);
    greenCone=thisFrame(:,:,2);
    blueCone=thisFrame(:,:,3);
    linesImage=zeros(size(thisFrame, 1), size(thisFrame, 2));
    if find(trackedFrames==itFrame)>0
        thisFrameRows=fullInfo(find(fullInfo(:,9)==trackedFrames(itFrame)),:);
        for itSegment=1:size(thisFrameRows, 1)
            linesImage=or(linesImage, func_Drawline(zeros(size(thisFrame, 1),...
                size(thisFrame, 2)), thisFrameRows(itSegment, 2),...
            thisFrameRows(itSegment, 1), thisFrameRows(itSegment, 7), ...
            thisFrameRows(itSegment, 6), 1));
            linesImage=addNumbersImage(linesImage, thisFrameRows(itSegment, 10), [7, 14],...
                [thisFrameRows(itSegment, 2), thisFrameRows(itSegment, 1)]);
        end
    end
       
    redCone(linesImage~=0)=linesImage(linesImage~=0)*255;
    greenCone(linesImage~=0)=linesImage(linesImage~=0)*255;
    
    rbgCone(:,:,1)=redCone;
    rbgCone(:,:,2)=greenCone;
    rbgCone(:,:,3)=blueCone;
    imwrite(rbgCone, [imagesFolder '\Segments\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'jpg')
end


%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'Segments\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');
close(itwaitbar)
itwaitbar = waitbar(1,'Finished with segment images');
pause(.1)
close(itwaitbar)
