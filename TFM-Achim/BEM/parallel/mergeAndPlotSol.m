function[fx fy x_out y_out sol_coef]=mergeAndPlotSol(path)
load([path,'/files/imputFile.mat']);

M=mergeSolutions(path,imputFile.numTasks);

[fx fy x_out y_out sol_coef]=calcSolFromFwdMap(M,imputFile.u,imputFile.forceMesh,imputFile.L);

quiver(x_out, y_out, fx, fy)

