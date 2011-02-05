function [constrForceField]=identifySpecialCells(constrForceField,imageFileList)
% This has to be executed before trackNetwork.
% This function only works on networks (clusters >=2).

if nargin <2 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.TIF;*.tif;*.jpg;*.png;*.*'}, ...
       'Select the first image (of e.g. myosin marker)');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'cell') && ~isempty(constrForceField{frame}.cell)
        toDoList=horzcat(toDoList,frame);
    end
end


for frame=toDoList
    I = double(imread(imageFileList{frame}));

    % first determine the average intensity within the cells
    numCells=length(constrForceField{frame}.cell);
    for cellID=1:numCells
        % integrate the intesity within the cell:
        iVal_raw(cellID) =sum(sum(I(constrForceField{frame}.cell{cellID}.innerMask)));
        cell_area(cellID)=sum(constrForceField{frame}.cell{cellID}.innerMask(:));
    end
    % integrate the intensity outside the cell:
    iVal_bg =sum(sum(I(~constrForceField{frame}.segmRes.mask)));
    area_bg =sum(~constrForceField{frame}.segmRes.mask(:));
    
    % 1) The most simple, and most of the time, the most effective one:
    %meanI=mean(I(:));
    
    % 2) The bg in the image (wo cells), this one doesn't perform well:
    %meanI=iVal_bg/area_bg;
    
    % 3) The bg in the cells:
    meanI=sum(iVal_raw)/sum(cell_area);
    
    % Calculate the intensity value:
    iVal(frame,1:numCells)=iVal_raw-meanI*cell_area;
    
%     meanI=mean(I(:));
%     for cellID=1:length(constrForceField{frame}.cell)
%         % integrate the intesity within the cell:
%         iVal_raw=sum(sum(I(constrForceField{frame}.cell{cellID}.innerMask)));
%         
%         % substract the average image intensity weighted with the cell
%         % area:
%         iVal(frame,cellID)=iVal_raw-meanI*sum(constrForceField{frame}.cell{cellID}.innerMask(:));
%     end
    clear cell_area iVal_raw     
end

% -------------------------------------------------------------------------
% Old code:
% -------------------------------------------------------------------------

% remove the ones that are eactly 0 (These are actually empty entries,
% since the field got extended when a new cell entered the cluster):
iVal(iVal==0)=NaN;


% Normalize the distribution, that all fall into [0,1]:
iVal=(iVal-min(iVal(:)))/(max(iVal(:))-min(iVal(:)));

% sort the intensity values, this should give (in the case of n cells)
% n-seperated distributions:
for i=1:length(iVal)
    iVal_sorted(i,:)=sort(iVal(i,:));
end

% show the sorted intensity distribution:
figure(1)
hist(iVal_sorted,linspace(0,1,101));
xlim([-0.05 1.05])

% ask the user to put in an appropriate threshold:
threshold  = input(['Set the threshold value (',num2str(length(toDoList)),'frames): ']);
markerChar = input('Is this nuclei marker for control (put 0) or myosin cells (put [1])');
if isempty(markerChar) || markerChar~=0
    markerChar=1;
end

myoType    = input('Which myosin type is it? type e.g. [myoIIA_hp93], myoIIA_hp94, myoIIB_hp103: ','s');
if isempty(myoType)
    myoType='myoIIA_hp93';
end

% store the values in the cell structure:
numSpec=0;
numNorm=0;
for frame=toDoList
    for cellID=1:length(constrForceField{frame}.cell)
        if (iVal(frame,cellID)>threshold && markerChar==1) || (iVal(frame,cellID)<threshold && markerChar==0)            
            constrForceField{frame}.cell{cellID}.stats.spec=1;
            constrForceField{frame}.cell{cellID}.stats.type=myoType;            
            numSpec=numSpec+1;
        else
            constrForceField{frame}.cell{cellID}.stats.spec=0;
            constrForceField{frame}.cell{cellID}.stats.type='ctrl';     
            numNorm=numNorm+1;
        end
    end    
end
display(['normal cells: ',num2str(numNorm),';   myosin cells: ',num2str(numSpec)]);

