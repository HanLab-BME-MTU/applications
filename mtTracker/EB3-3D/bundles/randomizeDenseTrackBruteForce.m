function [randTrack,firstSearch]=randomizeDenseTrackBruteForce(track,dynManifoldsCell,mappingDistance,maxRandomDist,varargin)
% randomize a track so that no tracks map in the manifolds

firstSearch=TracksHandle;

trackDiffCoeff=arrayfun(@(t) nanmean(sum([t.x(1)-t.x(2:end); t.y(1)-t.y(2:end); t.z(1)-t.z(2:end)].^2))/(6*t.lifetime) ,track);

manifoldCell=cell(1,length(dynManifoldsCell));
livingManifold=zeros(1,length(dynManifoldsCell));
fIdx=track.f(1);
for mIdx=1:length(manifoldCell)
  pMidx=find(dynManifoldsCell{mIdx}(1).f==fIdx);
  if(~isempty(pMidx));
      manifoldCell{mIdx}= [[dynManifoldsCell{mIdx}(1).x(pMidx); dynManifoldsCell{mIdx}(1).y(pMidx);dynManifoldsCell{mIdx}(1).z(pMidx)], ...
                           [dynManifoldsCell{mIdx}(2).x(pMidx); dynManifoldsCell{mIdx}(2).y(pMidx);dynManifoldsCell{mIdx}(2).z(pMidx)]];
      livingManifold(mIdx)=1;
  end
end
manifoldCell=manifoldCell(logical(livingManifold));

inManifold=true;
point=[track.x(1);track.y(1);track.z(1)];
searchPoint=[];
randInitCount=0;
while inManifold
  randDir=rand(3,1)-0.5;
  randPoint=point+maxRandomDist*randDir;
  inManifold=belongTo(randPoint,manifoldCell,mappingDistance);
  searchPoint= [searchPoint randPoint];
  randInitCount=randInitCount+1;
end

firstSearch=TracksHandle;
firstSearch.x=searchPoint(1,:);
firstSearch.y=searchPoint(2,:);
firstSearch.z=searchPoint(3,:);
firstSearch.startFrame=1;
firstSearch.endFrame=size(searchPoint,2);

randTrack=TracksHandle;
randTrack.x(1)=randPoint(1);
randTrack.y(1)=randPoint(2);
randTrack.z(1)=randPoint(3);
randTrack.startFrame=track.startFrame;
randTrack.endFrame=track.endFrame;

maxCount=1000;
randCount=0;
for pIdx=2:length(track.f)
  fIdx=track.f(pIdx);
  manifoldCell=cell(1,length(dynManifoldsCell));
  livingManifold=zeros(1,length(dynManifoldsCell));
  for mIdx=1:length(manifoldCell)
      pMidx=find(dynManifoldsCell{mIdx}(1).f==fIdx);
      if(~isempty(pMidx));
          manifoldCell{mIdx}= [[dynManifoldsCell{mIdx}(1).x(pMidx); dynManifoldsCell{mIdx}(1).y(pMidx);dynManifoldsCell{mIdx}(1).z(pMidx)], ...
              [dynManifoldsCell{mIdx}(2).x(pMidx); dynManifoldsCell{mIdx}(2).y(pMidx);dynManifoldsCell{mIdx}(2).z(pMidx)]];
          livingManifold(mIdx)=1;
      end
  end
    manifoldCell=manifoldCell(logical(livingManifold));
  inManifold=true;
  while (inManifold)&&(randCount<maxCount)
    R = normrnd(randPoint,sqrt(trackDiffCoeff));
    inManifold=belongTo(R,manifoldCell,mappingDistance);
    randCount=randCount+1;
  end
  randTrack.x(pIdx)=R(1);
  randTrack.y(pIdx)=R(2);
  randTrack.z(pIdx)=R(3);
end
disp(['Exploratory trials:: Init: ' num2str(randInitCount) ' tracks: ' num2str(randCount) '/' num2str(length(track.f)-1)]); 

%  randKin=randKinTracksISO(tIdx)
%   randKin.x=min(randKin.x,MD.imSize_(1));
%   randKin.x=max(randKin.x,1);
%   randKin.y=min(randKin.y,MD.imSize_(2));
%   randKin.y=max(randKin.y,1);
%   randKin.z=min(randKin.z,MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_);
%   randKin.z=max(randKin.z,1);


function mapped=belongTo(point,manifoldCell,mappingDistance)
  mIdx=1;
  mapped=false;
  while (~mapped)&&(mIdx<length(manifoldCell)+1)
      mapped=mapPointsTo1DManifold(point,manifoldCell{mIdx},mappingDistance);
      mIdx=mIdx+1;
  end
