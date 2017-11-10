function M=mergeSolutions(path,numSol)

load([path,'/files/sol_',num2str(1,'%03d')]);
M=sol;
clear('sol');
    
for j=2:numSol
    load([path,'/files/sol_',num2str(j,'%03d')]);
    M=M+sol;
    clear('sol');
end

save([path,'/files/completeSol'],'M');    