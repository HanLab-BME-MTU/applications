%% Data Preparation
meshPtsFwdSol=2^7;
xmax=2^7;
ymax=2^7;
zmax=2^4;
gridSpacing=1;
xmin=1;
ymin=1;
zmin=1;

[x_mat_u,y_mat_u,z_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax,zmin:gridSpacing:zmax);

E=8000;
forceType = 'groupForce';
opt=[];
h=5000;
v=0.5;
refine='true';
useSameSampling='false';
%% Original Force Generation
% posNA = [139.0000  267.0000 150.00
%          156.0000  232.0000 150.00
%          184.0000  212.0000 150.00
%          217.0000  200.0000 150.00
%          246.0000  195.0000 150.00
%          272.0000  199.0000 150.00
%          297.0000  211.0000 150.00
%          323.0000  231.0000 150.00
%          346.0000  256.0000 150.00];
%      
%      posFA = [   42.0000  323.0000 150.0000
%                  54.0000  319.0000 150.0000
%                  96.0000  298.0000 150.0000
%                 139.0000  288.0000 150.0000
%                 180.0000  280.0000 150.0000
%                 225.0000  276.0000 150.0000
%                 263.0000  273.0000 150.0000
%                 301.0000  275.0000 150.0000
%                 331.0000  281.0000 150.0000
%                 351.0000  279.0000 150.0000
%                 381.0000  290.0000 150.0000
%                 417.0000  303.0000 150.0000
%                 455.0000  317.0000 150.0000
%                  69.0000  359.0000 150.0000
%                  91.0000  349.0000 150.0000
%                 115.0000  337.0000 150.0000
%                 129.0000  338.0000 150.0000
%                 168.0000  327.0000 150.0000
%                 186.0000  321.0000 150.0000
%                 270.0000  318.0000 150.0000
%                 288.0000  320.0000 150.0000
%                 312.0000  323.0000 150.0000
%                 361.0000  330.0000 150.0000
%                 383.0000  337.0000 150.0000
%                 423.0000  352.0000 150.0000
%                 432.0000  363.0000 150.0000];
%             
%             posBigFA = [ 79    434  150
%                          145   409  150
%                          230   376  150
%                          306   384  150
%                          377   419  150
%                          438   447  150];
% force_x = zeros(size(x_mat_u));
% force_y = zeros(size(y_mat_u));
      force_xx = assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,139/2,257/2,15000,62000,5000,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,156/4,232/4,11000/2,65000/2,2000,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,184/4,212/4,6000/2,70000/2,2500,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,217/4,200/4,2000,72000,-2500,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,246/4,195/4,0,75000/2,2000,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,272/4,199/4,-3000/2,71000/2,15000,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,297/4,211/4,-7500,69500,10000,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,323/4,231/4,-10000/2,64000/2,-15000,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,346/4,256/4,-14000,60000,12500,400/72,500/72,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,42/4,323/4,68000,160000,25000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,54/4,319/4,60000,166000,-15000,500/108,2100/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,96/4,298/4,52000,172000,5000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,139/4,288/4,44000,178000,97500,350/108,2300/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,180/4,280/4,16000,184000,105000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,225/4,276/4,8000,190000,105000,300/108,2400/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,263/4,273/4,0,200000,65000,550/108,2300/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,301/4,275/4,-8000,190000,-65000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,331/4,281/4,-16000,184000,55000,500/108,2200/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,351/4,279/4,-24000,178000,75000,400/108,2000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,381/4,290/4,-32000,172000,-45000,600/108,2700/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,417/4,303/4,-40000,166000,-85000,500/108,2300/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,455/4,317/4,-48000,160000,-75000,450/108,2100/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,69/4,359/4,60000,250000,-30000,700/108,3000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,91/4,349/4,50000,260000,30000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,115/4,337/4,40000,270000,40000,700/108,3500/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,129/4,338/4,30000,280000,-50000,700/108,3000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,168/4,327/4,20000,290000,60000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,186/4,321/4,10000,300000,25000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,270/4,318/4,0,310000,30000,-150/108,3500/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,288/4,320/4,-10000,300000,-30000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,312/4,323/4,-20000,290000,20000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,361/4,330/4,-30000,280000,40000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,383/4,337/4,-40000,270000,35000,600/108,3500/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,423/4,352/4,-50000,260000,20000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(1,x_mat_u,y_mat_u,z_mat_u,432/4,363/4,-60000,250000,30000,600/108,2600/108,forceType);
            
      force_yy = assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,139/4,257/4,15000,62000,5000,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,156/4,232/4,11000/2,65000/2,2000,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,184/4,212/4,6000/2,70000/2,2500,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,217/4,200/4,2000,72000,-2500,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,246/4,195/4,0,75000/2,2000,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,272/4,199/4,-3000/2,71000/2,15000,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,297/4,211/4,-7500,69500,10000,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,323/4,231/4,-10000/2,64000/2,-15000,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,346/4,256/4,-14000,60000,12500,400/72,500/72,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,42/4,323/4,68000,160000,25000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,54/4,319/4,60000,166000,-15000,500/108,2100/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,96/4,298/4,52000,172000,5000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,139/4,288/4,44000,178000,97500,350/108,2300/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,180/4,280/4,16000,184000,105000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,225/4,276/4,8000,190000,105000,300/108,2400/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,263/4,273/4,0,200000,65000,550/108,2300/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,301/4,275/4,-8000,190000,-65000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,331/4,281/4,-16000,184000,55000,500/108,2200/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,351/4,279/4,-24000,178000,75000,400/108,2000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,381/4,290/4,-32000,172000,-45000,600/108,2700/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,417/4,303/4,-40000,166000,-85000,500/108,2300/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,455/4,317/4,-48000,160000,-75000,450/108,2100/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,69/4,359/4,60000,250000,-30000,700/108,3000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,91/4,349/4,50000,260000,30000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,115/4,337/4,40000,270000,40000,700/108,3500/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,129/4,338/4,30000,280000,-50000,700/108,3000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,168/4,327/4,20000,290000,60000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,186/4,321/4,10000,300000,25000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,270/4,318/4,0,310000,30000,-150/108,3500/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,288/4,320/4,-10000,300000,-30000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,312/4,323/4,-20000,290000,20000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,361/4,330/4,-30000,280000,40000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,383/4,337/4,-40000,270000,35000,600/108,3500/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,423/4,352/4,-50000,260000,20000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(2,x_mat_u,y_mat_u,z_mat_u,432/4,363/4,-60000,250000,30000,600/108,2600/108,forceType);
 force_x = zeros(128,128,16);
  force_x(:,:,8) = force_xx(:,:,8);
 force_y = zeros(128,128,16);
  force_y(:,:,8) = force_yy(:,:,8);          
            
% force_z = zeros(size(z_mat_u));
% force_z(:,:,8) = force_x(:,:,8).*force_y(:,:,8);
      force_zz = assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,139/4,257/4,15000,62000,5000,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,156/4,232/4,11000/2,65000/2,2000,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,184/4,212/4,6000/2,70000/2,2500,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,217/4,200/4,2000,72000,-2500,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,246/4,195/4,0,75000/2,2000,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,272/4,199/4,-3000/2,71000/2,15000,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,297/4,211/4,-7500,69500,10000,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,323/4,231/4,-10000/2,64000/2,-15000,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,346/4,256/4,-14000,60000,12500,400/72,500/72,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,42/4,323/4,68000,160000,25000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,54/4,319/4,60000,166000,-15000,500/108,2100/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,96/4,298/4,52000,172000,5000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,139/4,288/4,44000,178000,97500,350/108,2300/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,180/4,280/4,16000,184000,105000,500/108,2000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,225/4,276/4,8000,190000,105000,300/108,2400/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,263/4,273/4,0,200000,65000,550/108,2300/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,301/4,275/4,-8000,190000,-65000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,331/4,281/4,-16000,184000,55000,500/108,2200/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,351/4,279/4,-24000,178000,75000,400/108,2000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,381/4,290/4,-32000,172000,-45000,600/108,2700/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,417/4,303/4,-40000,166000,-85000,500/108,2300/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,455/4,317/4,-48000,160000,-75000,450/108,2100/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,69/4,359/4,60000,250000,-30000,700/108,3000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,91/4,349/4,50000,260000,30000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,115/4,337/4,40000,270000,40000,700/108,3500/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,129/4,338/4,30000,280000,-50000,700/108,3000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,168/4,327/4,20000,290000,60000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,186/4,321/4,10000,300000,25000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,270/4,318/4,0,310000,30000,-150/108,3500/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,288/4,320/4,-10000,300000,-30000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,312/4,323/4,-20000,290000,20000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,361/4,330/4,-30000,280000,40000,600/108,2600/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,383/4,337/4,-40000,270000,35000,600/108,3500/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,423/4,352/4,-50000,260000,20000,600/108,3000/108,forceType)+...
                assumedForceAniso3D(3,x_mat_u,y_mat_u,z_mat_u,432/4,363/4,-60000,250000,30000,600/108,2600/108,forceType);

%  f= -x_mat_u.^2 + y_mat_u.^2 +z_mat_u.^2;
%  [force_x8,force_y8,force_z8] = gradient(f);
%  force_x(:,:,8) = force_x8(:,:,8);
%  force_y(:,:,8) = force_y8(:,:,8);
 force_z = zeros(128,128,16);
   force_z(:,:,8) = force_zz(:,:,8);          
%    force_z(:,:,9) = force_zz(:,:,9);          
%    force_z(:,:,7) = force_zz(:,:,7);          

% 
% 
% 
% for i=1:((1/2)*(max(max(max(z_mat_u))))-1)
%  force_x(:,:,i) = zeros(size(x_mat_u(:,:,i)));
% end
% for i=((1/2)*(max(max(max(z_mat_u))))+1):max(max(max(z_mat_u)))
%  force_x(:,:,i) = zeros(size(x_mat_u(:,:,i)));
% end
% for i=1:((1/2)*(max(max(max(z_mat_u))))-1)
%  force_y(:,:,i) = zeros(size(x_mat_u(:,:,i)));
% end
% for i=((1/2)*(max(max(max(z_mat_u))))+1):max(max(max(z_mat_u)))
%  force_y(:,:,i) = zeros(size(x_mat_u(:,:,i)));
% end
% for i=1:((1/2)*(max(max(max(z_mat_u))))-1)
%  force_z(:,:,i) = zeros(size(x_mat_u(:,:,i)));
% end
% for i=((1/2)*(max(max(max(z_mat_u))))+1):max(max(max(z_mat_u)))
%  force_z(:,:,i) = zeros(size(x_mat_u(:,:,i)));
% end
%  force_x = 100000000.*force_x;
%  force_y = 100000000.*force_y;
%  force_z = 100000000.*force_z;
%% Displacement Field
[ux,uy,uz,x_grid,y_grid,z_grid] = fwdSolution3D(x_mat_u, y_mat_u,z_mat_u,E,xmin,xmax,ymin,ymax,zmin,zmax,force_x,force_y,force_z,'fft',opt,v);


%% Data Saving
save('uxuyuz.mat','ux','uy','uz','x_grid','y_grid','z_grid','meshPtsFwdSol')