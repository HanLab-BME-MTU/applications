function [Results]=LeaveOneOutRandResampTreeBagger(A,B,label)

Den = vertcat(A.den);
ind = sum(vertcat(A.comp),2)>12;
Den = Den(ind,:);
[Den,labels]=makeExtraFeatures(Den,label);

cellNumA = [];
  for i = 1:(numel(A))
	    cellNumA = vertcat(cellNumA,i*ones(size(A(i).den,1),1));
            i
end

	    cellNumA = cellNumA(ind);


	    cellNumB=[];

DenEGF = vertcat(B.den);
ind = sum(vertcat(B.comp),2)>12;
DenEGF = DenEGF(ind,:);
[DenEGF,labels]=makeExtraFeatures(DenEGF,label);
for i=1:numel(B)
	    cellNumB = vertcat(cellNumB,(i)*ones(size(B(i).den,1),1));
	    i
end

	    cellNumB = cellNumB(ind);
    

      Results = struct('VarImp',[],'varSC',[],'i',[]);

      compEGF=DenEGF;
      for j=1:numel(A)
       comp = Den;
       comp(cellNumA==j,:)=[];
       RandomResampleTreeBagger;	  
Results(j).VarImp=VarImp;
Results(j).varSc = varSc;
Results(j).i=j;
      end

      n = numel(A);   
comp=Den;
for k=1:numel(B)
     compEGF = DenEGF;
     compEGF(cellNumB==k,:)=[];
     j=k+n;
     RandomResampleTreeBagger;
     Results(k+n).VarImp=VarImp;
     Results(k+n).varSc = varSc;
     Results(k+n).i=k+n;
end

end

