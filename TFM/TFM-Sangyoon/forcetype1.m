imScale=1;
numPoints_u=50*imScale;       %(must be even number)
xmin =1;
xmax =50*imScale;
ymin =1;
ymax =50*imScale;

[x_mat_u, y_mat_u]=meshgrid(linspace(xmin,xmax,numPoints_u) , linspace(ymin,ymax,numPoints_u));
x_vec_u=reshape(x_mat_u,[],1);
y_vec_u=reshape(y_mat_u,[],1); % fine mesh for original force

fpos = [10*imScale,10*imScale;
        12*imScale,11*imScale;
        15*imScale,10*imScale;
        11*imScale,14*imScale; % 1st cluster
        
        15*imScale,26*imScale;
        17*imScale,31*imScale;
        28*imScale,27*imScale; % large adhesion
        21*imScale,20*imScale; % 2nd cluster
        
        27*imScale,39*imScale; 
        31*imScale,42*imScale; 
        37*imScale,43*imScale;
        31*imScale,36*imScale; %
        38*imScale,40*imScale;
        40*imScale,33*imScale; % 3rd cluster,  large adhesion

        33*imScale,9*imScale; 
        39*imScale,18*imScale; 
        43*imScale,24*imScale;
        38*imScale,12*imScale;
        42*imScale,17*imScale]; % 4th cluster
fvec = [1200,-100;
        1800,-160;
        1500,-130;
        1300,-110;
        
        1100,-380;
        1300,-480;
        850,-280; % large adhesion, small force
        400,-10;

        500,-600;
        450,-700;
        250,-900;
        800,-1200; 
        280,-800;
        700,-2000; % large adhesion, large force

        -450,500;
        -400,550;
        -350,400;
        -350,500;
        -150,600];

force_x1 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(1,1),fpos(1,2),fvec(1,1),fvec(1,2),forceType,'smallFA');
force_y1 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(1,1),fpos(1,2),fvec(1,1),fvec(1,2),forceType,'smallFA'); %group1
force_x2 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(2,1),fpos(2,2),fvec(2,1),fvec(2,2),forceType,'smallFA');
force_y2 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(2,1),fpos(2,2),fvec(2,1),fvec(2,2),forceType,'smallFA'); %group1
force_x3 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(3,1),fpos(3,2),fvec(3,1),fvec(3,2),forceType,'smallFA');
force_y3 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(3,1),fpos(3,2),fvec(3,1),fvec(3,2),forceType,'smallFA'); %group1
force_x4 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(4,1),fpos(4,2),fvec(4,1),fvec(4,2),forceType,'smallFA');
force_y4 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(4,1),fpos(4,2),fvec(4,1),fvec(4,2),forceType,'smallFA'); %group2
force_x5 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(5,1),fpos(5,2),fvec(5,1),fvec(5,2),forceType,'smallFA');
force_y5 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(5,1),fpos(5,2),fvec(5,1),fvec(5,2),forceType,'smallFA'); %group2
force_x6 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(6,1),fpos(6,2),fvec(6,1),fvec(6,2),forceType,'smallFA');
force_y6 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(6,1),fpos(6,2),fvec(6,1),fvec(6,2),forceType,'smallFA'); %group2
force_x7 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(7,1),fpos(7,2),fvec(7,1),fvec(7,2),forceType,'smallFA');
force_y7 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(7,1),fpos(7,2),fvec(7,1),fvec(7,2),forceType,'smallFA'); %group2, large one
force_x8 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(8,1),fpos(8,2),fvec(8,1),fvec(8,2),forceType,'smallFA');
force_y8 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(8,1),fpos(8,2),fvec(8,1),fvec(8,2),forceType,'smallFA'); %group2
force_x9 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(9,1),fpos(9,2),fvec(9,1),fvec(9,2),forceType,'smallFA');
force_y9 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(9,1),fpos(9,2),fvec(9,1),fvec(9,2),forceType,'smallFA'); %group2
force_x10 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(10,1),fpos(10,2),fvec(10,1),fvec(10,2),forceType,'smallFA');
force_y10 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(10,1),fpos(10,2),fvec(10,1),fvec(10,2),forceType,'smallFA'); %group1


force_x11 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(11,1),fpos(11,2),fvec(11,1),fvec(11,2),forceType,'smallFA');
force_y11 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(11,1),fpos(11,2),fvec(11,1),fvec(11,2),forceType,'smallFA'); %group1
force_x12 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(12,1),fpos(12,2),fvec(12,1),fvec(12,2),forceType,'smallFA');
force_y12 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(12,1),fpos(12,2),fvec(12,1),fvec(12,2),forceType,'smallFA'); %group1
force_x13 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(13,1),fpos(13,2),fvec(13,1),fvec(13,2),forceType,'smallFA');
force_y13 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(13,1),fpos(13,2),fvec(13,1),fvec(13,2),forceType,'smallFA'); %group1
force_x14 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(14,1),fpos(14,2),fvec(14,1),fvec(14,2),forceType,'smallFA');
force_y14 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(14,1),fpos(14,2),fvec(14,1),fvec(14,2),forceType,'smallFA'); %group2, large one
force_x15 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(15,1),fpos(15,2),fvec(15,1),fvec(15,2),forceType,'smallFA');
force_y15 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(15,1),fpos(15,2),fvec(15,1),fvec(15,2),forceType,'smallFA'); %group2
force_x16 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(16,1),fpos(16,2),fvec(16,1),fvec(16,2),forceType,'smallFA');
force_y16 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(16,1),fpos(16,2),fvec(16,1),fvec(16,2),forceType,'smallFA'); %group2
force_x17 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(17,1),fpos(17,2),fvec(17,1),fvec(17,2),forceType,'smallFA');
force_y17 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(17,1),fpos(17,2),fvec(17,1),fvec(17,2),forceType,'smallFA'); %group2
force_x18 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(18,1),fpos(18,2),fvec(18,1),fvec(18,2),forceType,'smallFA');
force_y18 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(18,1),fpos(18,2),fvec(18,1),fvec(18,2),forceType,'smallFA'); %group2
force_x19 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(19,1),fpos(19,2),fvec(19,1),fvec(19,2),forceType,'smallFA');
force_y19 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(19,1),fpos(19,2),fvec(19,1),fvec(19,2),forceType,'smallFA'); %group2

force_x = force_x1+force_x2+force_x3+force_x4+force_x5+force_x6+...
        force_x7+force_x8+force_x9+force_x10+force_x11+force_x12+...
        force_x13+force_x14+force_x15+force_x16+force_x17+force_x18+force_x19;
force_y = force_y1+force_y2+force_y3+force_y4+force_y5+force_y6+...
        force_y7+force_y8+force_y9+force_y10+force_y11+force_y12+...
        force_y13+force_y14+force_y15+force_y16+force_y17+force_y18+force_y19;

[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,[],[],[],[],...
    @(x,y) assumedForceShifted(1,x,y,10*imScale,10*imScale,1200,-100,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,12*imScale,11*imScale,1800,-160,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,15*imScale,10*imScale,1500,-130,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,11*imScale,14*imScale,1300,-110,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,15*imScale,26*imScale,1100,-380,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,17*imScale,31*imScale,1300,-480,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,28*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,40*imScale,33*imScale,700,-1000,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,33*imScale,9*imScale,-450,-800,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,39*imScale,18*imScale,-400,550,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,43*imScale,24*imScale,-350,400,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,38*imScale,12*imScale,-350,500,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,42*imScale,17*imScale,-150,600,forceType,'smallFA'),...
    @(x,y) assumedForceShifted(2,x,y,10*imScale,10*imScale,1200,-100,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,12*imScale,11*imScale,1800,-160,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,15*imScale,10*imScale,1500,-130,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,11*imScale,14*imScale,1300,-110,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,15*imScale,26*imScale,1100,-380,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,17*imScale,31*imScale,1300,-480,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,28*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,40*imScale,33*imScale,700,-1000,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,33*imScale,9*imScale,-450,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,39*imScale,18*imScale,-400,550,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,43*imScale,24*imScale,-350,400,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,12*imScale,-350,500,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,42*imScale,17*imScale,-150,600,forceType,'smallFA'),'fft',[],meshPtsFwdSol);
