function [displFieldNew]=prepDisplForBEM(displField,opt)
% prepare the displacements for BEM method such that:
% The points in the displ list are always at exactly the same position. In
% particular, the list has always the same length.
% If a measurement at a point is missing in one frame the user has the
% following options to fix this problem:
% INPUT
% opt: 'sortout': The point will be killed (sorted out) in all frames.
%      'linear' : This point is recovered by linear interpolation
%      'nearest': This point is recovered by nearest neighbor interpolation
% In case of 'linear' it might happen that points at the field boundary can
% not be interpolated and yield NaNs. Those are interpolated by nearest
% neighbor interpolations.
       


% all position entries in one column:
pos=vertcat(displField.pos);
maxX=max(pos(:,1));
maxY=max(pos(:,2));

% transform to linear indices
ind=sub2ind([maxY maxX],pos(:,2),pos(:,1));

% count multiple indices. Optimal would be that each index appears the same
% time (number of frames)
allPts=unique(ind); 
[counts]=histc(ind,allPts);

% The number of different positions ever used:
numPtsAll=length(allPts);

if numPtsAll~=length(counts)
    erorr('length of these two lists must agree!!!')
end

% kill all indices that don't appear in every frame:
maxOcc=max(counts);

if maxOcc>length(displField)
    erorr('The cut off level can (should) be the number of entries (frames) in the displField!!!')
end

% Find the positions with full/not-full occurance:
checkVec=(counts==maxOcc);

% the list of permanent points:
permPts=allPts(checkVec);

if strcmp(opt,'sortout')    
    numPtsNew=length(permPts);
    % take only those position in each frame that are found in the upper
    % list:
    [ypos, xpos]=ind2sub([maxY maxX],permPts);
       
    % find the corresponding displ-vecs:   
    for iframe=1:length(displField)
       matXcord=NaN*zeros(maxY,maxX);
       matYcord=NaN*zeros(maxY,maxX);
        
       indForCurrFrame=sub2ind([maxY maxX],displField(iframe).pos(:,2),displField(iframe).pos(:,1));
       matXcord(indForCurrFrame)=displField(iframe).vec(:,1);
       matYcord(indForCurrFrame)=displField(iframe).vec(:,2);
       
       % Now only pull the good values:
       displFieldNew(iframe).vec(:,1)=matXcord(permPts);
       displFieldNew(iframe).vec(:,2)=matYcord(permPts);
       
       % fill in the positions (this is kind of stupid, since they are the
       % same for all frames):
       displFieldNew(iframe).pos(:,1)=xpos;
       displFieldNew(iframe).pos(:,2)=ypos;
       
       %transfer also the parameter field:
       displFieldNew(iframe).par=displField(iframe).par;
       displFieldNew(iframe).par.prep4fastBEM=1;
       
       % show the number of killed positions:
       display(['Pts killed in frame ',num2str(iframe),': ',num2str(length(displField(iframe).pos(:,1))-numPtsNew)])
    end
else
    for iframe=1:length(displField)
       matXcord=NaN*zeros(maxY,maxX);
       matYcord=NaN*zeros(maxY,maxX);
       
       
       indForCurrFrame=sub2ind([maxY maxX],displField(iframe).pos(:,2),displField(iframe).pos(:,1));
       matXcord(indForCurrFrame)=displField(iframe).vec(:,1);
       matYcord(indForCurrFrame)=displField(iframe).vec(:,2);
       
       % find the points that are missing in this frame: 
       missPts=setdiff(allPts,indForCurrFrame);
       
       % transform back to x/y-coordinates:
       [ypos, xpos]=ind2sub([maxY maxX],missPts);
       
       % interpolate the displacement at these positions. Create the
       % interpolants:
       ux_intp=TriScatteredInterp(displField(iframe).pos(:,1),displField(iframe).pos(:,2),displField(iframe).vec(:,1),opt);
       uy_intp=TriScatteredInterp(displField(iframe).pos(:,1),displField(iframe).pos(:,2),displField(iframe).vec(:,2),opt);
       
       % Perform the interpolation:
       ux=ux_intp(xpos,ypos);
       uy=uy_intp(xpos,ypos);
       
       % Check that everything went well:
       bdPts=find(isnan(ux));
       if ~isempty(bdPts)
           % Create the nearest neighbor interpolants for the
           % boundary points:
           ux_intp=TriScatteredInterp(displField(iframe).pos(:,1),displField(iframe).pos(:,2),displField(iframe).vec(:,1),'nearest');
           uy_intp=TriScatteredInterp(displField(iframe).pos(:,1),displField(iframe).pos(:,2),displField(iframe).vec(:,2),'nearest');
           
           % Perform the interpolation:
           ux(bdPts)=ux_intp(xpos(bdPts),ypos(bdPts));
           uy(bdPts)=uy_intp(xpos(bdPts),ypos(bdPts));
       end
       
       % add those values to the matrix above:
       matXcord(missPts)=ux;
       matYcord(missPts)=uy;
       
       % Now read out the full list:
       displFieldNew(iframe).vec(:,1)=matXcord(allPts);
       displFieldNew(iframe).vec(:,2)=matYcord(allPts);       
       
       % Also fill in the positions:
       [yposAll, xposAll]=ind2sub([maxY maxX],allPts);
       
       displFieldNew(iframe).pos(:,1)=xposAll;
       displFieldNew(iframe).pos(:,2)=yposAll;
       
       %transfer also the parameter field:
       displFieldNew(iframe).par=displField(iframe).par;
       displFieldNew(iframe).par.prep4fastBEM=1;
       
       % show the number of interpolated positions:
       display(['Pts intp in frame ',num2str(iframe),': ',num2str(length(missPts))]);
    end
end
    

return;
% to compare the results:
opt='linear';
[displFieldNew]=prepDisplForBEM(displField,opt);

for iframe=1:length(displFieldNew)
    if strcmp(opt,'sortout')
        quiver(   displField(iframe).pos(:,2),   displField(iframe).pos(:,1),   displField(iframe).vec(:,2),   displField(iframe).vec(:,1),0,'r');
        hold on;
        quiver(displFieldNew(iframe).pos(:,2),displFieldNew(iframe).pos(:,1),displFieldNew(iframe).vec(:,2),displFieldNew(iframe).vec(:,1),0,'g');
        hold off;
    else
        quiver(displFieldNew(iframe).pos(:,2),displFieldNew(iframe).pos(:,1),displFieldNew(iframe).vec(:,2),displFieldNew(iframe).vec(:,1),0,'r');
        hold on;        
        quiver(   displField(iframe).pos(:,2),   displField(iframe).pos(:,1),   displField(iframe).vec(:,2),   displField(iframe).vec(:,1),0,'g');
        hold off;
    end
    input('enter to continue: ');
end