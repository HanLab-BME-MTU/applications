function info=analyzeSpeckles(imgSize,lMax,lMin)

% SYNOPSIS   info=analyzeSpeckles(imgSize,lMax,lMin)

% Setting vertex coordinates
dSet=[1,1;imgSize(1),1;imgSize(1),imgSize(2);1,imgSize(2)];

% Collect coordinates into matrices
[y x]=find(lMax);
pMax=[y x];
if isempty(pMax)
   info=[];
   return;
end
[y x]=find(lMin);
pMin=[y x];

% Attach dSet to lMin
pMin = cat(1,pMin,dSet);

% Delaunay triangulation
triMin=delaunay(pMin(:,1),pMin(:,2));

% Search triangles
triangles=tSearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));

% Store information
for i=1:size(triangles,1)
   if ~isnan(triangles(i))  % If NaN -> no triangle found
      info(1,1:2,i)=pMax(i,:);
      info(2:4,1,i)=pMin([triMin(triangles(i),:)],1);
      info(2:4,2,i)=pMin([triMin(triangles(i),:)],2);  
   else   % Actually remove locmax by setting Imax = Imin
      info(1,1:2,i)=pMax(i,:);
      info(2:4,1,i)=pMax(i,1);
      info(2:4,2,i)=pMax(i,2);
   end
end;
