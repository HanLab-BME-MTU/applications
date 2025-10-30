function force = assumedForceAniso3D(j,x,y,z,xshift,yshift,wx,wy,wz,d1,d2,forceType)
% force = assumedForceAniso2D(j,x,y,xshift,yshift,wx,wy,d1,d2,forceType) takes grid x and y and make anisotropic gaussian
% distributed force field of which sources are xshift and yshift.
% input     :   x           grid of x coordinates
%               y           grid of y coordinates
%               z           grid of z coordinates
%               xshift      x value of point source of force
%               yshift      y value of point source of force
%               wx          x component of force orientation
%               wy          y component of force orientation
%               wz          z component of force orientation
%               d1          diameter of adhesion in the direction of minor
%               axis
%               d2          diameter of adhesion in the direction of major
%               axis
%              
%               forceType   type of force
%               ('pointForce','groupForce', or 'smoothForce')
%               FAsize      'largeFA' or 'smallFA'
%               dx,dy       diameter of FA in pixel in directions of normal
%               and tangential to (wx, wy) assuming that the length of FA
%               is in parallel with the direction of force
%               
% output    :   force       x or y or z grid of force distribution (if j=1 or 2 or 3,
% respectively)
%              
% Sangyoon Han Jan 2013
% if strcmp(FAsize,'largeFA')
%     std = 4.5;
%     adh_r = 4.1; % adhesion diameter in pixel
% else
%     std = 1.5;
%     adh_r = 1.1;
% end

% orientation of force
theta = atan2(wy,wx);
% psi = atan2(wz,(wx^2+wy^2)^(1/2));
% if theta > pi/2

stdx = d1/2*1.1;
stdy = d2/2*1.1;
adh_rx = d1/2;
adh_ry = d2/2;

% z=z(1,1,:);
%  for i=1:length(z)
%      if abs(z(i)-8) < 0.5
%          disp('in the loop')

if j==1
    switch(forceType)
        case 'groupForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wx*aniso2DMask(xshift, yshift, adh_ry, adh_rx, theta, x, y).*force;
        case 'pointForce'
            force=wx*(x==xshift || y==yshift );
        case 'smoothForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wx*force;
    end
    elseif j==2
    switch(forceType)
        case 'groupForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wy*aniso2DMask(xshift, yshift, adh_ry, adh_rx, theta, x, y).*force;
        case 'pointForce'
            force=wy*(x==xshift || y==yshift);
        case 'smoothForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wy*force;
     
     end
    elseif j==3
       switch(forceType)
        case 'groupForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wz*aniso2DMask(xshift, yshift, adh_ry, adh_rx, theta, x, y).*force;
        case 'pointForce'
            force=wz*(x==xshift || y==yshift);
        case 'smoothForce'
            force = anisoGaussian2DMatrix(xshift, yshift, stdy, stdx, theta, x, y);
            % cutting out outside of area
            force = wz*force;
      end
else
    error('please input 1 or 2 or 3 for j');
end
return

%      end
%  end
end

