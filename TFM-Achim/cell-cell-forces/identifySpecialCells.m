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
    meanI=mean(I(:));
    for cellID=1:length(constrForceField{frame}.cell)
        % integrate the intesity within the cell:
        iVal_raw=sum(sum(I(constrForceField{frame}.cell{cellID}.innerMask)));
        
        % substract the average image intensity weighted with the cell
        % area:
        iVal(frame,cellID)=iVal_raw-meanI*sum(constrForceField{frame}.cell{cellID}.innerMask(:));
    end    
end

% -------------------------------------------------------------------------
% Old code:
% -------------------------------------------------------------------------

% Normaliz the distribution, that all fall into [0,1]:
iVal=iVal/max(iVal(:));

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
threshold=input('Set the threshold value: ');

% store the values in the cell structure:
for frame=toDoList
    for cellID=1:length(constrForceField{frame}.cell)
        if iVal(frame,cellID)>threshold
            constrForceField{frame}.cell{cellID}.stats.spec=1;
        else
            constrForceField{frame}.cell{cellID}.stats.spec=0;
        end
    end    
end

