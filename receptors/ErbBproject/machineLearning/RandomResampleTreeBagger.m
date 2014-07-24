% randomly samples type = 0 and type = 1 and runs TreeBagger 50 times
% uses 50 trees
% returns of 0.05 confidence levels by permutation
%
% 
%

VarImp = zeros(20,20);
varSc = zeros(50,20);
randVarImp = zeros(20,20);
randVarSc=varSc;

ParaOps =statset;
ParaOps.UseParallel = true;

%Take an equal number of type 0 and type 1


for i = 1:50
    n = fix(numel(comp(:,1))*0.66);
    ind = randperm(numel(comp(:,1)),n);
    indEGF = randperm(numel(compEGF(:,1)),n);
    
    type = vertcat(zeros(size(ind))',ones(size(indEGF))');
    measure = vertcat(comp(ind,:),compEGF(indEGF,:));
b = TreeBagger(50,measure,type,'OOBVarImp','On','MinLeaf',300,'Options',ParaOps);
    [t,idx] = sort(b.OOBPermutedVarDeltaError);
    idx2 = sub2ind([20,20],1:20,idx);
    VarImp(idx2)=VarImp(idx2)+1;
varSc(i,:)=b.OOBPermutedVarDeltaError;

    ind = randperm(numel(type));
b = TreeBagger(50,measure,type(ind),'OOBVarImp','On','MinLeaf',300,'Options',ParaOps);
    [t,idx] = sort(b.OOBPermutedVarDeltaError);
    idx2 = sub2ind([20,20],1:20,idx);
    randVarImp(idx2)=randVarImp(idx2)+1;
randVarSc(i,:)=b.OOBPermutedVarDeltaError;
i    
end
2
figure,subplot(1,2,1);
image(VarImp(20:-1:1,:)*3)
ylabel('Feature rank');
title(['Variable Importance Histogram (Starved v Stim) #',num2str(j)])
subplot(1,2,2)
image(randVarImp(20:-1:1,:)*3)
colormap('hot')
title('Variable Importance Histogram (Random Assignment)')
