function [newTrajMat,trajMat]=brTrajectoryModif(trajectoryDescription);
% Reshape the interpolation by recopying the values over the numb of frame
% with similar value



trajMat=cat(1,trajectoryDescription.individualStatistics.dataListGroup);
trajMat=trajMat(:,[1:5 7]);
%newTrajMat=trajMat;

trajMat=trajMat(find(trajMat(:,3)==1 |trajMat(:,3)==2),:);

testDiff=trajMat(:,2)-trajMat(:,1);

newMatIndex=1;

for i=1:length(testDiff)
    if testDiff(i)==1  
        newTrajMat(newMatIndex,:)=trajMat(i,:);
    elseif testDiff(i)>1   
        for j=newMatIndex:newMatIndex+testDiff(i)-1
            newTrajMat(j,1)=trajMat(i,1)+j-newMatIndex;
            newTrajMat(j,2)=newTrajMat(j,1)+1;
            newTrajMat(j,3:end)=trajMat(i,3:end);
            newTrajMat(j,6)=trajMat(i,6)/testDiff(i);
            
        end 
        
    end
    newMatIndex=newMatIndex+testDiff(i);
end

        