function []=drawSheetEdge(imageFileListPhaseInput,toDoListInput)
% This is a largely simplified version of the drawSheetEdge. The old
% version (which is unstable) can be found in the deprecated folder!

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


if nargin<2 || isempty(toDoListInput)
    toDoListInput=1:length(imageFileListPhase);
end

replyReuse=0;
allFiles=dir;
for fileId=1:length(allFiles)
    if strcmp(allFiles(fileId).name,'xDetectEdgeHand.mat')
        replyReuse=input('Found a file with hand drawn edges. Use it? Yes=1, No=0. [Yes]:');
    end
end


if isempty(replyReuse) || replyReuse==1
    % Overwrite the automatically detected values by the ones found in:
    load('xDetectEdgeHand.mat');
end


for frame=toDoListInput
    [maxRows,maxCols]=size(sheetMask(frame).mat);
    markerPix=sheetMask(frame).mat(dPix(frame)+1,dPix(frame)+1);
    
    drawByHand=0;
    oldBnD =sheetBnD(frame).pos;
    oldEdge=sheetEdge(frame).pos;
    oldMask=sheetMask(frame).mat;
    
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
        sheetMask(frame).mat =oldMask;
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
        sheetMask(frame).mat =sheetMaskMat;
        sheetBnD(frame).pos  =sheetBnDpos;
        sheetEdge(frame).pos =fliplr(sheetEdgePos.coord);
        
        %         figure(2)
        %         imagesc(sheetMask(frame).mat)
        %         hold on
        %         plot(sheetBnD(frame).pos(:,2),sheetBnD(frame).pos(:,1),'k')
        %         plot(sheetEdge(frame).pos(:,2),sheetEdge(frame).pos(:,1),'.r')
        %         hold off
    end
    toDoList=union(toDoList,frame);
    save('xDetectEdgeHand.mat', 'sheetMask', 'sheetEdge', 'sheetBnD','toDoList','-v7.3');
end
display('Done! You should now re-run cellMig_part_2_statistics for this data set!')