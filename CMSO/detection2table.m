function T=detection2table(detections,timeInterval,channel)
% P. Roudot 2016
varName={'C','t','X','Y','Z','Int','dX','dY','dZ','dInt'};
tableCell=cell(1,length(detections));
for fIdx=1:length(detections)
    if(isfield(detections(fIdx),'zCoord'))
        Z=detections(fIdx).zCoord(:,1);
        dZ=detections(fIdx).zCoord(:,2);
    else
        Z=zeros(size((detections(fIdx).xCoord(:,1))));
        dZ=Z;
    end
    detArray=[channel*ones(size(Z)), timeInterval*ones(size(Z)), detections(fIdx).xCoord(:,1), detections(fIdx).yCoord(:,1),Z,detections(fIdx).int(:,1), detections(fIdx).xCoord(:,2), detections(fIdx).yCoord(:,2),dZ(:),detections(fIdx).int(:,2)];
    tableCell{fIdx}=array2table(detArray,'variableNames',varName);
end
T=table();
for fIdx=1:length(detections)
    T=[T;tableCell{fIdx}];
end