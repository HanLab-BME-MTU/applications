function sol=performTask(path,j)
load([path,'/files/imputFile.mat']);

sol=calcFwdMap(imputFile.x, imputFile.y, imputFile.forceMesh, imputFile.E,  imputFile.task(j).span);

save([path,'/files/sol_',num2str(j, '%03d')],'sol');