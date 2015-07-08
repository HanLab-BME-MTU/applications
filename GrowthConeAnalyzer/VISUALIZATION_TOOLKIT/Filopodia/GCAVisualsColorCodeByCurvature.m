function [ output_args ] = plotByCurvature( filoInfo )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [ny,nx] = size(img); 
    
  filoCurvsAll =   abs(vertcat(filoInfo(:).Ext_FiloCurvature));
  xy = vertcat(filoInfo(:).Ext_coordsXY_SplineFit); 
% create distance mapper 
cMapLength=128; cMap=jet(cMapLength);
            mapper=linspace(min(filoCurvsAll),max(filoCurvsAll),cMapLength)'; 

 D=createDistanceMatrix(filoCurvsAll,mapper);
            [sD,idxCMap]=sort(abs(D),2);
setFigure(nx,ny,'on')
imshow(-img,[]) 
hold on 

              for k=1:cMapLength
%                 idxCand = EPreFilt(idxCMap(:,1) == k,2);
%                 idxSeed = EPreFilt(idxCMap(:,1)==k,1);
%                 for iEdge = 1:length(idxCand) % some can have the same color
                    scatter(xy(idxCMap(:,1)==k,1),xy(idxCMap(:,1) == k,2),10,...
                      cMap(k,:),'filled');
                %end
             end

end

