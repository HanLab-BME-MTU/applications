function rotMat=mkRotMat(e_1,e_2,e_3)
%creates rotation matrices that project e1,e2,e3 onto 100,010,001
%
%SYNOPSIS rotMat=mkRotMat(e_1,e_2,e_3)
%
%INPUT lists of unit vectors e1-3 spanning 3D spaces
%
%OUTPUT rotMat: rotational matrix to turn a point in e1-3 to 100,010,001
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nVectors=size(e_1,1);

rotMat=zeros(3,3,nVectors);

for i=1:nVectors
    
    if ~all(e_1(i,:)==0);
        %rotMat =[e1,e2,e3]' since all three are unit vectors
        rotMat(:,:,i)=[e_1(i,:);e_2(i,:);e_3(i,:)];
    end
end