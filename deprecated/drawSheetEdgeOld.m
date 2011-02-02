function []=drawSheetEdge(imageFileListPhaseInput,toDoListInput)

try
    load('xDetectEdge.mat');
    load('xDetectNuclei.mat');  % for dPix
    load('xParameters.mat')
    % This overloads the old toDolist
    load('xResults.mat')
catch exception
    display('Couldnt find all essential files, please browse to folder')
    return;
end
  
if nargin<1 || isempty(imageFileListPhaseInput)
    imageFileListPhase=[pwd,filesep,'Phase'];
    try
        imageFileListPhase=getFileListFromFolder(imageFileListPhase);
    catch exception2
        [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
            'Select First phase Image');
        
        if ~ischar(filename) || ~ischar(pathname)
            return;
        end
        
        imageFileListPhase = getFileStackNames([pathname filesep filename]);
    end
elseif isdir(imageFileListPhaseInput)
    imageFileListPhase=getFileListFromFolder(imageFileListPhaseInput);
else
    isValid = 1;
    for frame = 1:numel(imageFileListPhaseInput)
        isValid = isValid && exist(imageFileListPhaseInput{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
    imageFileListPhase=imageFileListPhaseInput;
end
firstRun=1;
doMore=0;

if nargin<2 || isempty(toDoListInput)
    lastGood=toDoList(end);
    toDoListRest=(lastGood+1:length(sheetMask));
    if isempty(toDoListRest)
        toDoListRest=(lastGood+1:length(imageFileListPhase));
        doMore  =1;
        firstRun=0;
    end
else
    display('To-do-list extended to length of sheetMask.')
    toDoListRest=[toDoListInput , toDoListInput(end)+1:length(sheetMask)];
    lastGood=toDoListInput(1)-1;        
end

sheetMaskOld=sheetMask;
sheetBnDOld=sheetBnD;
sheetEdgeOld=sheetEdge;

clear sheetMask sheetEdge sheetBnD toDoList

maxHandFrames=0;
try
    load('xDetectEdgeHand.mat');
    replyReuse=input('Found a file with hand drawn edges. Use it? Yes=1, No=0. [Yes]:');
    if isempty(replyReuse) || replyReuse==1
        % Overwrite the automatically detected values by the ones found in
        % the hand drawn file:
        maxHandFrames=length(sheetMask);
        sheetMaskOld(1:maxHandFrames)=sheetMask;
        sheetBnDOld(1:maxHandFrames) =sheetBnD;
        sheetEdgeOld(1:maxHandFrames)=sheetEdge;

%!!!    This line is bad. The whole thing should be moved up.
        toDoListRest=maxHandFrames+1:length(imageFileListPhase);
        
        clear sheetMask sheetEdge sheetBnD toDoList;
    else
        clear sheetMask sheetEdge sheetBnD toDoList;
    end
end

% sheetMask(1:lastGood)=sheetMaskOld(1:lastGood);
% sheetBnD(1:lastGood) =sheetBnDOld(1:lastGood);
% sheetEdge(1:lastGood)=sheetEdgeOld(1:lastGood);

sheetMask=sheetMaskOld;
sheetBnD =sheetBnDOld;
sheetEdge=sheetEdgeOld;


drawByHand(1:lastGood)=0;


while firstRun==1 || doMore==1
    if firstRun==1 && doMore==0
        for frame=toDoListRest
            if length(sheetMaskOld)<=frame
                display('Restart the program with: "drawSheetEdge;" to draw remaining frames.')
                return;
            end
            [maxRows,maxCols]=size(sheetMaskOld(frame).mat);                
            markerPix=sheetMaskOld(frame).mat(dPix(frame)+1,dPix(frame)+1);
            
            drawByHand=0;
            oldBnD =sheetBnDOld(frame).pos;
            oldEdge=sheetEdgeOld(frame).pos;
            oldMask=sheetMaskOld(frame).mat;
            
            currentImage=double(imread(imageFileListPhase{frame}));
            
            rectangle=zeros([maxRows,maxCols]);
            rectangle((1+dPix(frame)):(maxRows-dPix(frame)),(1+dPix(frame)):(maxCols-dPix(frame)))=1;
            B = bwboundaries(rectangle,4);
            rectCurve=B{1};
            
            figure(1)
            imagesc(currentImage)
            colormap('gray')
            hold on
            plot(oldBnD(:,2),oldBnD(:,1),'k');
            plot(oldEdge(:,2),oldEdge(:,1),'b');
            plot(rectCurve(:,2),rectCurve(:,1),'k');
            hold off
            title(['Old boundary frame: ',num2str(frame)]);
            axis equal
            xlim([1 maxCols])
            ylim([1 maxRows])
            
            replyDraw=input('Draw this sheet edge by hand? (No=enter; 1= take old edge as start; 2= start from scratch;)? [No]: ');
            if isempty(replyDraw) || replyDraw==0
                drawByHand(frame)=0;
                sheetMask(frame).mat=oldMask;
                sheetBnD(frame).pos  =oldBnD;
                sheetEdge(frame).pos =oldEdge;
            else
                drawByHand(frame)=replyDraw;
            end
            
            skipPix=25;
            
            if drawByHand(frame)>0
                curveBnD=[];
                while size(curveBnD,1)<2
                    display('Draw the edge now!')
                    if drawByHand(frame)==1
                        polygonObject = impoly(gca,fliplr(oldEdge(1:skipPix:end,:)),'Closed',false);
                        replyDummy=input('Press ENTER to continue:...');
                    elseif drawByHand(frame)==2
                        polygonObject = impoly(gca,'Closed',false);
                    else
                        display('Something went wrong!');
                    end
                    curveBnD       = round(getPosition(polygonObject));
                end
                curveBnD(curveBnD(:,2)<1,2)=1;
                curveBnD(curveBnD(:,2)>maxRows,2)=maxRows;
                % devide the mask into two parts:
                [sheetMask1, sheetMask2, sheetBnD1, sheetBnD2, sheetEdgePos, cell1_extMask, ~, ~, ~]=intersecMaskPolygon(fliplr(rectCurve),curveBnD);
                % extend the masks to the size of the image:
                sheetMask1(maxRows,maxCols)=0;
                sheetMask2(maxRows,maxCols)=0;
                
                % pick the right mask:
                if sheetMask1(dPix(frame)+1,dPix(frame)+1)==markerPix
                    sheetMaskMat=sheetMask1;
                    sheetBnDpos=fliplr(sheetBnD1);
                elseif sheetMask2(dPix(frame)+1,dPix(frame)+1)==markerPix
                    sheetMaskMat=sheetMask2;
                    sheetBnDpos=fliplr(sheetBnD2);
                else
                    display('Something went wrong');
                    return;
                end
                sheetMask(frame).mat=sheetMaskMat;
                sheetBnD(frame).pos  =sheetBnDpos;
                sheetEdge(frame).pos =fliplr(sheetEdgePos.coord);
                
                %         figure(2)
                %         imagesc(sheetMask(frame).mat)
                %         hold on
                %         plot(sheetBnD(frame).pos(:,2),sheetBnD(frame).pos(:,1),'k')
                %         plot(sheetEdge(frame).pos(:,2),sheetEdge(frame).pos(:,1),'.r')
                %         hold off
            end
            toDoList=1:max(frame,maxHandFrames);
            save('xDetectEdgeHand.mat', 'sheetMask', 'sheetEdge', 'sheetBnD','toDoList','-v7.3');
        end        
    elseif firstRun==0 || doMore==1
        for frame=toDoListRest
            % here we have to draw by hand, no sheet is provided!

            currentImage=double(imread(imageFileListPhase{frame}));
            [maxRows,maxCols]=size(currentImage);

            rectangle=zeros([maxRows,maxCols]);
            rectangle((1+dPix(frame)):(maxRows-dPix(frame)),(1+dPix(frame)):(maxCols-dPix(frame)))=1;
            B = bwboundaries(rectangle,4);
            rectCurve=B{1};

            figure(1)
            imagesc(currentImage)
            colormap('gray')
            hold on
            plot(rectCurve(:,2),rectCurve(:,1),'k');
            hold off
            title(['Old boundary frame: ',num2str(frame)]);
            axis equal
            xlim([1 maxCols])
            ylim([1 maxRows])

            drawByHand(frame)=1;
            curveBnD=[];
            while size(curveBnD,1)<2
                display('Draw the edge now!')
                polygonObject = impoly(gca,'Closed',false);
                curveBnD       = round(getPosition(polygonObject));
            end
            curveBnD(curveBnD(:,2)<1,2)=1;
            curveBnD(curveBnD(:,2)>maxRows,2)=maxRows;
            % devide the mask into two parts:
            [sheetMask1, sheetMask2, sheetBnD1, sheetBnD2, sheetEdgePos, cell1_extMask, ~, ~, ~]=intersecMaskPolygon(fliplr(rectCurve),curveBnD);
            % extend the masks to the size of the image:
            sheetMask1(maxRows,maxCols)=0;
            sheetMask2(maxRows,maxCols)=0;
            
            % pick the right mask (by largest overlap with previous mask):
            checkMat1 = sheetMask1 & sheetMask(frame-1).mat;
            checkMat2 = sheetMask2 & sheetMask(frame-1).mat;
            
            if sum(checkMat1(:))>sum(checkMat2(:))
                sheetMaskMat=sheetMask1;
                sheetBnDpos=fliplr(sheetBnD1);
            else
                sheetMaskMat=sheetMask2;
                sheetBnDpos=fliplr(sheetBnD2);
            end
            sheetMask(frame).mat=sheetMaskMat;
            sheetBnD(frame).pos  =sheetBnDpos;
            sheetEdge(frame).pos =fliplr(sheetEdgePos.coord);
            
            toDoList=1:max(frame,maxHandFrames);
            save('xDetectEdgeHand.mat', 'sheetMask', 'sheetEdge', 'sheetBnD','toDoList','-v7.3');
            
%             figure(2)
%             imagesc(sheetMask(frame).mat)
%             hold on
%             plot(sheetBnD(frame).pos(:,2),sheetBnD(frame).pos(:,1),'k')
%             plot(sheetEdge(frame).pos(:,2),sheetEdge(frame).pos(:,1),'.r')
%             hold off
        end
    end
    toDoList=1:max(frame,maxHandFrames);
    save('xDetectEdgeHand.mat', 'sheetMask', 'sheetEdge', 'sheetBnD','toDoList','-v7.3');
    
    firstRun=0;
    if toDoListRest(end)<length(imageFileListPhase)
        doMore=input('Do you think there might be more frames that could be analyzed 0=no; 1=yes: [no]');
        if isempty(doMore)
            doMore=0;
        end
    end
end
display('Done! You should now re-run cellMig_part_2_statistics for this data set!')