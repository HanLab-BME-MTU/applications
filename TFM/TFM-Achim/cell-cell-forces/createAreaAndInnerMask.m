function [constrForceField]=createAreaAndInnerMask(constrForceField)
% The following lines are needed to find the inner masks:
toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'cell') && ~isempty(constrForceField{frame}.cell)
        toDoList=horzcat(toDoList,frame);
    end
end

for frame=toDoList
    [rowsROI,colsROI]=size(constrForceField{frame}.segmRes.mask);    
    for cellID=1:length(constrForceField{frame}.cell)
        cellMask=constrForceField{frame}.cell{cellID}.mask;
        cellMask(rowsROI,colsROI)=0;        
        
        constrForceField{frame}.cell{cellID}.innerMask = constrForceField{frame}.segmRes.mask &  cellMask;
        constrForceField{frame}.cell{cellID}.cellArea  = sum(constrForceField{frame}.cell{cellID}.innerMask(:))*(constrForceField{frame}.par.pixSize_mu^2);
    end
end