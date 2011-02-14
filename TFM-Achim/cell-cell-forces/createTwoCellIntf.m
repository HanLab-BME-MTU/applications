function constrForceField=createTwoCellIntf(constrForceField,frame,doPlot)
% This function breaks drawn interfaces apart such that each interface is
% only connecting one cell pair: a two cell interface. Note that drawn
% interfaces might go along the interfaces of different cell pairs. The code
% is so bulky since we test for several special cases which arise from the
% fact that all interfaces are pixel-paths. The code makes sure that
% two-cell-interfaces never intersect, there are no stupid dump pieces
% originating at the intersection of pixel paths, interfaces never end in
% the limbo but reach out till either a vertex in the cluster or till the
% dilated curve surrounding the cluster is reached.

i=frame;
% first delete all existing twoCellIntf:
constrForceField{i}.twoCellIntf=[];

numItf=1;
for j=1:length(constrForceField{i}.interface)
    % intialization:
    m=1;
    ptsCurrIntf=length(constrForceField{i}.interface{j}.conct(:,1));
    pattern=constrForceField{i}.interface{j}.conct(m,:);
    
    % run through the whole length of the interface:
    while m<=ptsCurrIntf
        numPos=1;    
        % store the values for the first point, the starting vertex:
        twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
        twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);

        % determine which cells are connected through this interface. The
        % following procedure ensures that we get correct results if the
        % interface is composed of only two highly connected vertices:
        neigh1=find(constrForceField{i}.interface{j}.conct(m  ,:));
        neigh2=find(constrForceField{i}.interface{j}.conct(m+1,:));

        % The interface links the intersection of these two sets:
        twoCellIntf{numItf}.link=intersect(neigh1,neigh2);

        % now we can go to the next point of the interface:
        m=m+1;
        numPos=numPos+1;

        % while the connectivity of the points is <3 or we are not yet at the
        % end of this interface we go on:
        
        pattern=constrForceField{i}.interface{j}.conct(m,:);    
        while sum(pattern)<3 && m<=ptsCurrIntf            
            twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
            twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);
            m=m+1;
            numPos=numPos+1;
            if m<=ptsCurrIntf
                pattern=constrForceField{i}.interface{j}.conct(m,:);
            end
        end
        % check why the loop has terminated:
        if sum(pattern)>2 && m<ptsCurrIntf
            % we arrived at a higher order vertex within the interface:
            twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
            twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);
            
            % Don't increase m, the next interface should start at the same
            % vertex!
            numItf=numItf+1;
        elseif sum(pattern)>2 && m==ptsCurrIntf
            % we arrived at the end of the interface which is a higher 
            % order vertex. The position m=ptsCurrIntf has still to be
            % treated:
            twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
            twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);
            
            % Now the interface has been run through completely, to 
            % terminate the loop:
            m=m+1;
            numItf=numItf+1;
        else
            % we arrived at the end of the interface which is a vertex of
            % connectivity<3. The interface has been run through
            % completely, to terminate the loop:
            m=m+1;
            numItf=numItf+1;
            % m is already at its maximum
        end        
    end
end

%**************************************************************************
% check if there are close by vertices that should be linked:
%**************************************************************************
anzItf=length(twoCellIntf);
allItf=1:anzItf;
for numItf=allItf
    % this is the start point:
    endPt(numItf,:)       =twoCellIntf{numItf}.pos(1,:);
    % this is the end point:
    endPt(anzItf+numItf,:)=twoCellIntf{numItf}.pos(end,:);
end
allEndPts=1:2*anzItf;
connMat = -ones(2*anzItf,2*anzItf);
for idEndPt=allEndPts
    % exclude the own interfaces:
    for idCmpPt=setdiff(allEndPts,[idEndPt anzItf+idEndPt])
        [~,diff]=compPts(endPt(idEndPt,:),endPt(idCmpPt,:));
        connMat(idEndPt,idCmpPt)=norm(diff);
    end
end

connMatOrg=connMat;
connMat(connMat>2)=-1;

% Now run through all lines and find the rows that have 1 but no 0
for row=1:2*anzItf
    badVert=find(connMat(row,:)==1);
    if ~isempty(badVert) && sum(connMat(row,:)==0)==0
        % check if the neighboring vertices are all the same. This should
        % always be true. The problem we are resolving here only occurs
        % where three interfaces come together (in particular, there are no
        % vertices with 4 interfaces).
        checkSum=0;
        for nghV1=badVert
            for nghV2=setxor(badVert,nghV1)
                if ~compPts(endPt(nghV1,:),endPt(nghV2,:))              
                    checkSum=checkSum+1;
                end
            end
        end
        % Check if the interface can be shortened, if not extend it
        found=0;
        if ~checkSum
            if row<=anzItf
                currIntf=twoCellIntf{row}.pos;
                intfL=length(currIntf(:,1));
                k=1;
                while ~found && k<min(7,intfL)
                    if compPts(currIntf(k,:),endPt(badVert(1),:))
                        twoCellIntf{row}.pos(1:k-1,:)=[];
                        display(['Resolved bad vertex: frame: ',num2str(frame),'; intf: ',num2str(row),'; at k: ',num2str(k),'!']);
                        found=1;
                        % remove the resolved vertex from the matrix
                        for v=1:length(badVert)
                            connMat(row,badVert(v))=0;
                            % as well as the symmetric entry:
                            connMat(badVert(v),row)=0;
                        end                            
                    else
                        k=k+1;
                    end
                end
            else
                currIntf=twoCellIntf{row-anzItf}.pos;
                intfL=length(currIntf(:,1));
                k=intfL;
                while ~found && k>max([intfL-7,1])
                    if compPts(currIntf(k,:),endPt(badVert(1),:))
                        twoCellIntf{row-anzItf}.pos(k+1:end,:)=[];
                        display(['Resolved bad vertex: frame: ',num2str(frame),'; intf: ',num2str(row-anzItf),'; at k: ',num2str(k),'!']);
                        found=1;
                        % remove the resolved vertex from the matrix
                        for v=1:length(badVert)
                            connMat(row,badVert(v))=0;
                            % as well as the symmetric entry:
                            connMat(badVert(v),row)=0;
                        end
                    else
                        k=k-1;
                    end
                end
            end            
        else
            error('Dont know what to do. The checksum should always =0...')
        end
        if ~found
            error('Couldnt resolve bad vertex')        
        end
    end
end

% the suspicion margin:
dMargin=3;
if sum(sum(connMat==1))>0
    error('Not all bad vertices have been cleaned');
elseif connMatOrg>1 & connMatOrg<dMargin
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    display('!!!Some suspicious vertices have been left over!!!')
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
end


if nargin>2 && doPlot==1
    marker=['r','b','m','c','g','y','k'];
    figure(111)
    for p=1:length(twoCellIntf)
        plot(twoCellIntf{p}.pos(:,1),twoCellIntf{p}.pos(:,2),marker(mod(p,7)+1))
        hold on    
    end
    plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
    set(gca,'YDir','reverse')
    hold off
end

%**************************************************************************
% now check for singular interfaces (exist only outside the cell cluster)
%**************************************************************************
singList=[];
normList=[];
for numItf=1:length(twoCellIntf)
    intf     = twoCellIntf{numItf}.pos;
    mask     = constrForceField{frame}.segmRes.mask; % this is the inner mask that has 0 at holes!
    indIntf  = sub2ind(size(mask),intf(:,2), intf(:,1));
    checkVec = mask(indIntf);
    
    if sum(checkVec)==0
        twoCellIntf{numItf}.class='singular';
        singList=horzcat(singList,numItf);
    else
        twoCellIntf{numItf}.class='normal';
        normList=horzcat(normList,numItf);
    end    
end

% Now check which of the singular interfaces ends on the dilated mask.
% These are sorted into the singStrtList. The other ones start and end in
% the limbo (these are actually critical piece which can be linked
% incorrectly very easily)! The latter are sorted into singIntrList.
singStrtList=[];
singIntrList=[];
curveDilated=constrForceField{frame}.segmRes.curveDilated;
maxX=max(curveDilated(:,1));
maxY=max(curveDilated(:,2));
indCurv  = sub2ind([maxY maxX],curveDilated(:,2), curveDilated(:,1));
dilMask  = zeros(maxY,maxX);
dilMask(indCurv)  = 1;

for singIntf=singList
    intf     = twoCellIntf{singIntf}.pos;    
    indIntf  = sub2ind([maxY maxX],intf(:,2), intf(:,1)); 
    checkVec = dilMask(indIntf);
    
    if sum(checkVec)>0
        checkVec;
        twoCellIntf{singIntf}.class='singularStrt';
        singStrtList=horzcat(singStrtList,singIntf);
    else
        checkVec;
        twoCellIntf{singIntf}.class='singularIntr';
        singIntrList=horzcat(singIntrList,singIntf);
    end
end

% start with one out of singStrtList and connect it to all its potential
% neighbors. Check the status of the neighbor. If the neighbor is a normal
% interface, fuse it and change status to extended. If the neighbor is a
% singular edge, fuse it and put it back into the singStrtList. After
% running thorugh all interfaces, take the current one out of the
% singStrtList. Run through the singStrtList till empty!

searchList=1:length(twoCellIntf);
while ~isempty(singStrtList)
    strtIntf=singStrtList(1);
    searchList=setxor(searchList,strtIntf);
    
    foundList=[];
    for numItf=searchList
        
        strPt1=twoCellIntf{strtIntf}.pos(1  ,:);
        endPt1=twoCellIntf{strtIntf}.pos(end,:);
        
        strPt2=twoCellIntf{numItf}.pos(1  ,:);
        endPt2=twoCellIntf{numItf}.pos(end,:);
        
        if     compPts(strPt1,strPt2)
            % associate this interface with the neighboring one and
            % avoid duplicating the common point:
            twoCellIntf{numItf}.pos  = vertcat(flipud(twoCellIntf{strtIntf}.pos),twoCellIntf{numItf}.pos(2:end,:));
            if strcmp(twoCellIntf{numItf}.class,'normal') || strcmp(twoCellIntf{numItf}.class,'extended')
                twoCellIntf{numItf}.class='extended';
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat([NaN NaN],twoCellIntf{numItf}.pos);                
            elseif strcmp(twoCellIntf{numItf}.class,'singularIntr')
                twoCellIntf{numItf}.class='singularStrt';
                singStrtList=horzcat(singStrtList,numItf);
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat([NaN NaN],twoCellIntf{numItf}.pos);    
            else
                error('Something went wrong')
            end
            foundList = horzcat(foundList,strtIntf);
            
        elseif compPts(strPt1,endPt2)
            % associate this interface with the neighboring one and
            % avoid duplicating the common point:
            twoCellIntf{numItf}.pos  = vertcat(twoCellIntf{numItf}.pos(1:end-1,:),twoCellIntf{strtIntf}.pos);
            if strcmp(twoCellIntf{numItf}.class,'normal') || strcmp(twoCellIntf{numItf}.class,'extended')
                twoCellIntf{numItf}.class='extended';
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat(twoCellIntf{numItf}.pos,[NaN NaN]);    
            elseif strcmp(twoCellIntf{numItf}.class,'singularIntr')
                twoCellIntf{numItf}.class='singularStrt';
                singStrtList=horzcat(singStrtList,numItf);
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat(twoCellIntf{numItf}.pos,[NaN NaN]); 
            else
                error('Something went wrong')
            end
            foundList = horzcat(foundList,strtIntf);
            
        elseif compPts(endPt1,strPt2)
            % associate this interface with the neighboring one and
            % avoid duplicating the common point:
            twoCellIntf{numItf}.pos  = vertcat(twoCellIntf{strtIntf}.pos,twoCellIntf{numItf}.pos(2:end,:));
            if strcmp(twoCellIntf{numItf}.class,'normal') || strcmp(twoCellIntf{numItf}.class,'extended')
                twoCellIntf{numItf}.class='extended';
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat([NaN NaN],twoCellIntf{numItf}.pos); 
            elseif strcmp(twoCellIntf{numItf}.class,'singularIntr')
                twoCellIntf{numItf}.class='singularStrt';
                singStrtList=horzcat(singStrtList,numItf);
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat([NaN NaN],twoCellIntf{numItf}.pos); 
            else
                error('Something went wrong')
            end
            foundList = horzcat(foundList,strtIntf);
            
        elseif compPts(endPt1,endPt2)
            % associate this interface with the neighboring one and
            % avoid duplicating the common point:
            twoCellIntf{numItf}.pos  = vertcat(twoCellIntf{numItf}.pos(1:end-1,:),flipud(twoCellIntf{strtIntf}.pos));
            if strcmp(twoCellIntf{numItf}.class,'normal') || strcmp(twoCellIntf{numItf}.class,'extended')
                twoCellIntf{numItf}.class='extended';
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat(twoCellIntf{numItf}.pos,[NaN NaN]); 
            elseif strcmp(twoCellIntf{numItf}.class,'singularIntr')
                twoCellIntf{numItf}.class='singularStrt';
                singStrtList=horzcat(singStrtList,numItf);
                % cap the actual start point of that interface to avoid
                % double linking:
                twoCellIntf{numItf}.pos=vertcat(twoCellIntf{numItf}.pos,[NaN NaN]); 
            else
                error('Something went wrong')
            end
            foundList = horzcat(foundList,strtIntf);
            
        end
    end
    if ~isempty(foundList)
        singStrtList = setxor(singStrtList,foundList);
        twoCellIntf{foundList(1)}=[];
        twoCellIntf{foundList(1)}.class=[];
        
    else
        error('A starting singular interface could not be connected!?')
    end
end

% Now condense the field to avoid empty interfaces:
indx=0;
checkSum=0;
for numItf=1:length(twoCellIntf)
    if strcmp(twoCellIntf{numItf}.class,'extended') || strcmp(twoCellIntf{numItf}.class,'normal')
        indx=indx+1;
        checkVec=isnan(twoCellIntf{numItf}.pos(:,1));
        twoCellIntf{numItf}.pos(checkVec,:)=[];
        twoCellIntfCond{indx}=twoCellIntf{numItf};
    elseif ~isempty(twoCellIntf{numItf}.class)
        checkSum=checkSum+1;
    end
end

if checkSum~=0
    error('One singular interface could not be resolved!?')
end

twoCellIntf=twoCellIntfCond;

% now find interfaces that are only of length 3 and have a partner to match
% up with. This happens in the following situation:
%
%      s   s/* + + +
%+ + + s/*   *
%      o
%      o
% the 3 pixel long curves that are indivated by s and * pose the problem.
% They are created by dissecting the + and the o curves:
candidates=[];
for itfId=1:length(twoCellIntf)
    currL=size(twoCellIntf{itfId}.pos,1);
    if currL==3
        display(['Frame: ',num2str(frame),' itfId: ',num2str(itfId), ' might be bad edge, since length is only: ',num2str(currL)]);
        candidates=union(candidates,itfId);
    end
end

allItfIds=1:length(twoCellIntf);
allCandidates=candidates;
if length(candidates)>1
    % Then there might be a pair that we can resolve:
    % check if the start or endpoints match:
    while ~isempty(candidates)
        itfId1=candidates(1);
        candidates(1)=[];
        strPt1=twoCellIntf{itfId1}.pos(1,:);
        endPt1=twoCellIntf{itfId1}.pos(end,:);
        for itfId2=candidates
            strPt2=twoCellIntf{itfId2}.pos(1,:);
            endPt2=twoCellIntf{itfId2}.pos(end,:);
            if         (compPts(strPt1,strPt2) && compPts(endPt1,endPt2))...
                    || (compPts(strPt1,endPt2) && compPts(endPt1,strPt2))
                % Then the two small interfaces are essentially the
                % same (except the center point).
                % find the three fold vertex and the interface that
                % belongs to it. Then add one fragment from above to
                % this longer interface and kill the short
                % fragments.
                listStr=[];
                listEnd=[];
                for goodItfId=setdiff(allItfIds,allCandidates)
                    goodStrPt=twoCellIntf{goodItfId}.pos(1,:);
                    goodEndPt=twoCellIntf{goodItfId}.pos(end,:);
                    if compPts(strPt1,goodStrPt) || compPts(strPt1,goodEndPt)
                        listStr=horzcat(listStr,goodItfId);
                    elseif compPts(endPt1,goodStrPt) || compPts(endPt1,goodEndPt)
                        listEnd=horzcat(listEnd,goodItfId);
                    end
                end
                
                if length(listStr)==1 && length(listEnd)==2
                    % then it is a three fold vertx (since itfId1/2 also
                    % end there together with a good/long interface). The
                    % other vertex should be deg=4;
                    goodItfId=listStr;
                    goodStrPt=twoCellIntf{goodItfId}.pos(1,:);
                    goodEndPt=twoCellIntf{goodItfId}.pos(end,:);
                    if compPts(strPt1,goodStrPt)
                        input('Check this case!') % if the interfaces match and that no points are doubled!
                        twoCellIntf{goodItfId}.pos=vertcat(flipud(twoCellIntf{itfId1}.pos),twoCellIntf{goodItfId}.pos(2:end,:));
                        candidates=setdiff(candidates,[itfId1 itfId2]);
                        twoCellIntf{itfId1}=[];
                        twoCellIntf{itfId2}=[];
                        display(['Resolved 3-pix interfaces: ',num2str([itfId1 itfId2])]);
                    elseif compPts(strPt1,goodEndPt)
                        % ok!
                        twoCellIntf{goodItfId}.pos=vertcat(twoCellIntf{goodItfId}.pos(1:(end-1),:),twoCellIntf{itfId1}.pos);
                        candidates=setdiff(candidates,[itfId1 itfId2]);
                        twoCellIntf{itfId1}=[];
                        twoCellIntf{itfId2}=[];
                        display(['Resolved 3-pix interfaces: ',num2str([itfId1 itfId2])]);
                    else
                        input('This case should never happen, problems might occure later on.')
                    end
                    
                elseif length(listEnd)==1 && length(listStr)==2
                    goodItfId=listEnd;
                    goodStrPt=twoCellIntf{goodItfId}.pos(1,:);
                    goodEndPt=twoCellIntf{goodItfId}.pos(end,:);
                    if compPts(endPt1,goodStrPt)
                        % ok!
                        twoCellIntf{goodItfId}.pos=vertcat(twoCellIntf{itfId1}.pos,twoCellIntf{goodItfId}.pos(2:end,:));
                        candidates=setdiff(candidates,[itfId1 itfId2]);
                        twoCellIntf{itfId1}=[];
                        twoCellIntf{itfId2}=[];
                        display(['Resolved 3-pix interfaces: ',num2str([itfId1 itfId2])]);
                    elseif compPts(endPt1,goodEndPt)
                        input('Check this case!') % if the interfaces match and that no points are doubled!
                        twoCellIntf{goodItfId}.pos=vertcat(twoCellIntf{goodItfId}.pos(1:(end-1),:),flipud(twoCellIntf{itfId1}.pos));
                        candidates=setdiff(candidates,[itfId1 itfId2]);
                        twoCellIntf{itfId1}=[];
                        twoCellIntf{itfId2}=[];
                        display(['Resolved 3-pix interfaces: ',num2str([itfId1 itfId2])]);
                    end
                    
                else
                    input('This case should never happen, problems might occure later on.')
                end
            end
            
        end
    end
end
if ~isempty(candidates)
    input('There are still edges of length 3!')
end

% Now condense the field to avoid empty interfaces:
indx=1;
clear twoCellIntfCond;
for numItf=1:length(twoCellIntf)
    if ~isempty(twoCellIntf{numItf})        
        twoCellIntfCond{indx}=twoCellIntf{numItf};
        indx=indx+1;
    end    
end


marker=['r','b','m','c','g','y','k'];
figure(112)
for p=1:length(twoCellIntfCond)
    plot(twoCellIntfCond{p}.pos(:,1),twoCellIntfCond{p}.pos(:,2),marker(mod(p,7)+1))
    hold on
end
plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
set(gca,'YDir','reverse')
hold off

constrForceField{i}.twoCellIntf=twoCellIntfCond;

