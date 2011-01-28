function [constrForceFieldUpdated]=updateConstrForceField(constrForceField,forceField,method)
% updates the constrForceField give a new forceField. This function is
% useful e.g. when a forceField has been recalculated with a different
% regularization parameter or the displacement field has been reanalyzed
% with a different template size.

if nargin < 3 || isempty(method)
    method='pixIntp';
end

% check if the tow force fields are of the same length:
length_old=length(constrForceField);
length_new=length(forceField);

if length_old~=length_new
    numToDo=min(length_old,length_new);
    display('The lengths of old and new froce field dont match!')
else
    numToDo=length_old;
end

toDoList=[];
for frame=1:length(constrForceField)
    if ~isempty(constrForceField{frame}) && frame<=numToDo && isfield(constrForceField{frame},'cell');
        % only then, it is a good frame:
        toDoList=horzcat(toDoList,frame);
    end
end

for frame=toDoList
    ImgSize=size(constrForceField{frame}.segmRes.maskDilated);
    % displayText=[' ',num2str(frame),' of ',num2str(length(toDoList))];
    % progressText(frame/length(toDoList),displayText);
    % The fields: .network
    %             .network_tracked
    %need to be updated later on using the appropriate functions.
    
    % Fill in the unchanged fields:
    constrForceFieldUpdated{frame}.segmRes  =constrForceField{frame}.segmRes;
    constrForceFieldUpdated{frame}.interface=constrForceField{frame}.interface;
    if isfield(constrForceFieldUpdated{frame},'twoCellIntf');
        constrForceFieldUpdated{frame}.twoCellIntf=constrForceField{frame}.twoCellIntf;
    end
    
    
    %**********************************************************************
    % .par: fill in the new values                                        *
    %**********************************************************************
    constrForceFieldUpdated{frame}.par=forceField(frame).par;
    
    %**********************************************************************
    % .roi: fill in the new values                                        *
    %**********************************************************************
    fpos=forceField(frame).pos;
    fvec=forceField(frame).vec;
    if sum(sum(abs(fpos-round(fpos))))==0
        [inpos,invec]=findVectorFieldInMask(fpos,fvec,constrForceFieldUpdated{frame}.segmRes.maskDilated);
        
        constrForceFieldUpdated{frame}.roi.pos=inpos;
        constrForceFieldUpdated{frame}.roi.vec=invec;
    else
        display(['Use inpolygon for frame: ',num2str(frame)]);
        checkVector = inpolygon(fpos(:,1),fpos(:,2),constrForceFieldUpdated{frame}.segmRes.curveDilated(:,1),constrForceFieldUpdated{frame}.segmRes.curveDilated(:,2));
        
        constrForceFieldUpdated{frame}.roi.pos=fpos(checkVector,:);
        constrForceFieldUpdated{frame}.roi.vec=fvec(checkVector,:);
    end
    
    %**********************************************************************
    % .cell: fill in the new values                                        *
    %**********************************************************************
    for cellID=1:length(constrForceField{frame}.cell)
        
        % The fields: .uFwdSol
        %             .iuExp
        % need to be updated later on using the appropriate functions.
        
        % Fill in the unchanged fields:
        constrForceFieldUpdated{frame}.cell{cellID}.mask     =constrForceField{frame}.cell{cellID}.mask;
        constrForceFieldUpdated{frame}.cell{cellID}.extMask  =constrForceField{frame}.cell{cellID}.extMask;
        constrForceFieldUpdated{frame}.cell{cellID}.innerMask=constrForceField{frame}.cell{cellID}.innerMask;
        constrForceFieldUpdated{frame}.cell{cellID}.center   =constrForceField{frame}.cell{cellID}.center;
        constrForceFieldUpdated{frame}.cell{cellID}.boundary =constrForceField{frame}.cell{cellID}.boundary;
        constrForceFieldUpdated{frame}.cell{cellID}.interface=constrForceField{frame}.cell{cellID}.interface;
        constrForceFieldUpdated{frame}.cell{cellID}.innerMask=constrForceField{frame}.cell{cellID}.innerMask;
        constrForceFieldUpdated{frame}.cell{cellID}.cellArea =constrForceField{frame}.cell{cellID}.cellArea;
        
        %******************************************************************
        % .pos, vec and stats: fill in the new values                     *
        %******************************************************************
        
        
        % if the forces are defined at integer position we simply have
        % to find those that are located in the mask. Else one has to
        % use the slow inpolygon:
        ROIpos=constrForceFieldUpdated{frame}.roi.pos;
        ROIvec=constrForceFieldUpdated{frame}.roi.vec;
        if sum(sum(abs(ROIpos-round(ROIpos))))==0
            [inpos,invec]=findVectorFieldInMask(ROIpos,ROIvec,constrForceFieldUpdated{frame}.cell{cellID}.mask);
            
            constrForceFieldUpdated{frame}.cell{cellID}.pos=inpos;
            constrForceFieldUpdated{frame}.cell{cellID}.vec=invec;
        else
            checkVector = inpolygon(ROIpos(:,1),ROIpos(:,2),constrForceFieldUpdated{frame}.cell{cellID}.boundary(:,1),constrForceFieldUpdated{frame}.cell{cellID}.boundary(:,2));
            
            constrForceFieldUpdated{frame}.cell{cellID}.pos=ROIpos(checkVector,:);
            constrForceFieldUpdated{frame}.cell{cellID}.vec=ROIvec(checkVector,:);
        end
        
        % The stats field has many fields that need to be updated later on.
        constrForceFieldUpdated{frame}.cell{cellID}.stats.resForce.pos = constrForceFieldUpdated{frame}.cell{cellID}.center;
        
        
        %!!!!   % Here we multiply each stress with its support (gridsize^2*pixSize^2*10^(-12)) to
        % get the actual force in N. In the future, the support might
        % be triangles and then this becomes more complicated in case of no
        % interpolation! In case the force field is interpolated,
        % everything is fine also for irregular meshes. The conversion of
        % Pa to nN is performed by the function integrateForceField:
        pixSize_mu = constrForceFieldUpdated{frame}.par.pixSize_mu;
        bwMask = constrForceFieldUpdated{frame}.cell{cellID}.mask;
        if strcmp(method,'noIntp')
            gridSpacing= constrForceFieldUpdated{frame}.par.gridSpacing;
            [sumForceVec,method,~,~]=integrateForceField(forceField(frame).pos,forceField(frame).vec,bwMask,pixSize_mu,gridSpacing);
        else
            [sumForceVec,method,~,~]=integrateForceField(forceField(frame).pos,forceField(frame).vec,bwMask,pixSize_mu);
        end
        
        constrForceFieldUpdated{frame}.cell{cellID}.stats.resForce.vec = - sumForceVec;
        constrForceFieldUpdated{frame}.cell{cellID}.stats.resForce.mag =   sqrt(sum((constrForceFieldUpdated{frame}.cell{cellID}.stats.resForce.vec).^2));
        constrForceFieldUpdated{frame}.cell{cellID}.stats.method = method;
        
        display(['frame ',num2str(frame),', cell ',num2str(cellID),':']);
        display(['old force: ',num2str(constrForceField{frame}.cell{cellID}.stats.resForce.vec)]);
        display(['new force: ',num2str(constrForceFieldUpdated{frame}.cell{cellID}.stats.resForce.vec)]);
        display('----------------------------');
    end
    % Sum up the force over the whole cell cluster to get an error for each
    % force measurement:
    pixSize_mu = constrForceFieldUpdated{frame}.par.pixSize_mu;
    bwMask = constrForceFieldUpdated{frame}.segmRes.maskDilated;
    if strcmp(method,'noIntp')
        gridSpacing= constrForceFieldUpdated{frame}.par.gridSpacing;
        [errorSumForce,method,~,~]=integrateForceField(forceField(frame).pos,forceField(frame).vec,bwMask,pixSize_mu,gridSpacing);
    else
        [errorSumForce,method,~,~]=integrateForceField(forceField(frame).pos,forceField(frame).vec,bwMask,pixSize_mu);
    end
    
    constrForceFieldUpdated{frame}.errorSumForce.vec    = errorSumForce;
    constrForceFieldUpdated{frame}.errorSumForce.mag    = sqrt(sum((errorSumForce).^2));
    constrForceFieldUpdated{frame}.errorSumForce.method = method;
    
    display(['frame ',num2str(frame),', cell ',num2str(cellID),':']);
    display(['old error vector:        ',num2str(constrForceField{frame}.errorSumForce.vec)]);
    display(['new error vector:        ',num2str(constrForceFieldUpdated{frame}.errorSumForce.vec)]);
    display('----------------------------');
    
    % To check this (It works perfect!):
    showCheckPlot=0;
    if showCheckPlot==1
        figure(frame)
        marker=['r','b','m','c','g','y'];
        startPos=[0,0];
        for cellIndex=1:length(constrForceField{frame}.cell)
            currForce=constrForceFieldUpdated{frame}.cell{cellIndex}.stats.resForce.vec;
            plot([startPos(1) startPos(1)+currForce(1)],[startPos(2) startPos(2)+currForce(2)],marker(mod(cellIndex,6)+1));
            hold on;
            startPos=startPos+currForce;
        end
        plot([startPos(1) startPos(1)+errorSumForce(1)],[startPos(2) startPos(2)+errorSumForce(2)],'k','LineWidth',3);
        hold off
    end
    
    
end
display('All done!')
display('You will have to redo all steps after TFM_part_4_cutOutForceField!')

