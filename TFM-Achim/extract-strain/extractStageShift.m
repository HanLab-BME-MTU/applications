% Extract stage shift using the mpm matrix, stored in /track/mpm.mat
load('mpm.mat')
myM=M(:,:,24);
myM(~isfinite(1./myM(:,1)),:)=[];
myM(~isfinite(1./myM(:,3)),:)=[];
zeroLessM=myM;
mean(zeroLessM(:,1)-zeroLessM(:,3))
mean(zeroLessM(:,2)-zeroLessM(:,4))

