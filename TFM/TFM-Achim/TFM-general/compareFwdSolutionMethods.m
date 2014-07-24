function [iux_fast, iuy_fast, ux1_conv, uy1_conv]=compareFwdSolutionMethods(basisClass,x_vec,y_vec,E,meshPtsFwdSol)

% To make sure that the range over which the solution is calculated,
% take the double of the initial x and y ranges:
xmin=min(x_vec(:)); xmax=max(x_vec(:));
ymin=min(y_vec(:)); ymax=max(y_vec(:));
dx=xmax-xmin;
dy=ymax-ymin;

xrange=[-dx dx]';
yrange=[-dy dy]';


% Integration bounds:
xbd_min=min(basisClass.neighPos(:,1));
xbd_max=max(basisClass.neighPos(:,1));
ybd_min=min(basisClass.neighPos(:,2));
ybd_max=max(basisClass.neighPos(:,2));

% first basis solution, method = fft:
display('Calculate fast solution: ')
tic
[ux1_fast, uy1_fast, x_grid1_fast, y_grid1_fast]=fwdSolution(xrange,yrange,E,xbd_min,xbd_max,ybd_min,ybd_max,basisClass.basisFunc(1).f_intp_x,basisClass.basisFunc(1).f_intp_y,'fft','noIntp',meshPtsFwdSol);

% Then the interpolants of the first function are:
iux_fast = interp2(x_grid1_fast, y_grid1_fast, ux1_fast, x_vec, y_vec,'*cubic');
iuy_fast = interp2(x_grid1_fast, y_grid1_fast, uy1_fast, x_vec, y_vec,'*cubic');
toc

% Second basis solution, method = conv:
display('Calculate direct convolution: ')
[ux1_conv, uy1_conv]=fwdSolution(x_vec,y_vec,E,xbd_min,xbd_max,ybd_min,ybd_max,basisClass.basisFunc(1).f_intp_x,basisClass.basisFunc(1).f_intp_y,'conv');

rel_error_x=2*abs(ux1_conv-iux_fast)./abs(ux1_conv+iux_fast);
rel_error_y=2*abs(uy1_conv-iuy_fast)./abs(uy1_conv+iuy_fast);

mean_rel_error_x=mean(rel_error_x(:));
mean_rel_error_y=mean(rel_error_y(:));

figure(10000)
surf(x_vec,y_vec,rel_error_x)

figure(20000)
surf(x_vec,y_vec,rel_error_y)

display(['mean rel. error in x: ',num2str(mean_rel_error_x)]);
display(['max. rel. error in x: ',num2str(max(rel_error_x(:)))]);
display(['mean rel. error in y: ',num2str(mean_rel_error_y)]);
display(['max. rel. error in y: ',num2str(max(rel_error_y(:)))]);

return;

% To test the upper function use these lines:
xvec1=linspace(-500,500,2);
xvec2=linspace(-8,8,3);
xvec=sort(horzcat(xvec1, xvec2));

yvec1=linspace(-500,500,2);
yvec2=linspace(-8,8,3);
yvec=sort(horzcat(yvec1, yvec2));

[xvec yvec]=meshgrid(xvec,yvec);
E=10000;
meshPtsFwdSol=2^11;
[iux_fast, iuy_fast, ux1_conv, uy1_conv]=compareFwdSolutionMethods(basisClass,xvec,yvec,E,meshPtsFwdSol);

