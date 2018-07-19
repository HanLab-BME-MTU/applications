function constrForceField=createConctField(constrForceField,frame)
% find the neighboring domains to each point along the interfaces:
for j=1:length(constrForceField{frame}.interface)
    % The next line is actually not needed. The interfaces should have no
    % doubled points!
    constrForceField{frame}.interface{j}.pos=removeDoublePoints(constrForceField{frame}.interface{j}.pos);    
    
    currIntrf=constrForceField{frame}.interface{j}.pos;
    % initialize the connection field:
    constrForceField{frame}.interface{j}.conct=[];
    conct=[];
    for k=1:length(constrForceField{frame}.cell)
        currMask =constrForceField{frame}.cell{k}.extMask;
        conct(:,k)=cutIntfWithMask(currMask,currIntrf);
    end
    % clean up the connectivity matrix
    [conctCleaned, cleaned]=cleanUpConct(conct);
    if cleaned
        display(['Cleaned frame: ',num2str(frame),' interface: ',num2str(j)])
        % display(num2str([conct conctCleaned]));
    end
    constrForceField{frame}.interface{j}.conct=conctCleaned;
end

function [mat,didClean]=cleanUpConct(mat)
[rows,~]=size(mat);

checkVec= sum(abs(mat(2:rows,:)-mat(1:rows-1,:)),2)==0 & sum(mat(2:end,:),2)>2;
idBadPts=find(checkVec);
idBadPts=idBadPts(:)';
if ~isempty(idBadPts)
    didClean=true;
    
    % the connectivity values at these positions have to be set to the 2 value pattern!
    for badPt=idBadPts
        % get the 2 value pattern:
        currPos=badPt-1;
        found=false;
        shift=false;
        while ~found && currPos>=1
            if sum(mat(currPos,:),2)==2
                found=true;
                twoPattern=mat(currPos,:);
            else
                currPos=currPos-1;
            end
        end
        
        if ~found
            % Then we have reached the start of the interface without
            % finding a proper value. In this case we have to search in the
            % other direction
            
            %reinitialize the current position: 
            currPos=badPt+1;
            while ~found && currPos<=rows
                if sum(mat(currPos,:),2)==2
                    found=true;
                    shift=true;
                    twoPattern=mat(currPos,:);
                else
                    currPos=currPos+1;
                end
            end
        end
        
        if found
            if shift
                % if we had to go upwards we have to clean the following
                % position:
                mat(badPt+1,:)=twoPattern;
            else
                mat(badPt,:)=twoPattern;
            end
        else
            display('Interfaces are strange there might be an error!')
        end
        
    end
                
else
    didClean=false;
end
    