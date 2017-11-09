function [Results]=FeatureSelection_twoWays(data,type)

% assumes two classses only
  numF = size(data,2);

%use multilinear regression classifier

classResult = zeros(10,numF);
scramResult = [];

nZero = sum(type==0);
listZ = find(type==0);
nOne = sum(type==1);
listO = find(type);

n = fix(.66*min(nZero,nOne));
typeT = [zeros(n,1);ones(n,1)];

totalError = [];

for j=1:numF
	classlist = 1:numF;
classlist(j)=[];
for i=1:10
	indZ = randperm(nZero,n);
        indO = randperm(nOne,n);
        ind = [listZ(indZ),listO(indO)];
        dataT = data(ind,:);
        dataR = data;
        dataR(ind,:)=[];
        typeR = type;
        typeR(ind)=[];

%with all features
       % mdl = fitglm(dataT,typeT);
mdl = ClassificationKNN.fit(dataT,typeT,'NumNeighbors',20,'NSMethod','kdtree','Distance','euclidean');
        [pred,sc] = mdl.predict(dataR);
%pred(pred<0.5)=0;
%pred(pred>=0.5)=1;
        totalError = [totalError,sum(pred==typeR)/numel(pred)];

%scramble each feature
scram = zeros(1,numF);
parfor k=1:numF
	ind = randperm(size(dataR,1));
        temp = dataR;
        temp(:,k)=dataR(ind,k);
        pred = mdl.predict(temp);
        %pred(pred<0.5)=0;
        %pred(pred>=0.5)=1;
        scram(k)= (sum(pred==typeR)/numel(pred))-totalError(end);
end
scramResult = vertcat(scramResult,scram);       
display(['Left out #',num2str(j),' for the ',num2str(i),' time']);

%LeaveOneOut
        %mdl = fitglm(dataT(:,classlist),typeT);
        %pred = mdl.predict(dataR(:,classlist));
         %pred(pred<0.5)=0;
         %pred(pred>=0.5)=1;
mdl = ClassificationKNN.fit(dataT(:,classlist),typeT,'NumNeighbors',20,'NSMethod','kdtree','Distance','euclidean');
        [pred,sc] = mdl.predict(dataR(:,classlist));
classResult(i,j)=(sum(pred==typeR)/numel(pred))-totalError(end);
end
end


Results = struct('LeaveOneFeatureOut',classResult,'LeaveNoFeaturesOut',totalError,'ScrambleOneFeat',scramResult);

end
