function c = cCoefficientLagrangePlaneStress(E, nu, loc, state)
%cCoefficientLagrangePlaneStress Calculate c-coefficient for nonlinear plane stress
% Calculate the c-coefficient for a geometrically nonlinear Lagrangian formulation 
% of plane stress elasticity. The strain measure is the Green-Lagrange strain
% tensor. The stress is the second Piola-Kirchoff stress tensor. The material
% is assumed to be isotropic with linear behavior (Hooke's law applies).
% 
% E  - Young's modulus of the linear isotropic material
% nu - Poisson's ratio for the material
% p  - matrix of point (node) locations
% t  - element connectivity matrix
% u  - current displacement vector

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    31-Jan-2014 09:50:09

% Copyright 2014-2015 The MathWorks, Inc.


ux = reshape(state.ux,2,[]);
uy = reshape(state.uy,2,[]);



dudx=ux(1,:); 
dvdx=ux(2,:);
dudy=uy(1,:);  
dvdy=uy(2,:);
 
% if(~isempty(u))
%   [ux,uy] = pdegrad(p,t,u);
%   dudx=ux(1,:); dudy=uy(1,:); dvdx=ux(2,:); dvdy=uy(2,:);
% else
%   dudx = zeros(1, size(t,2)); dudy=dudx; dvdx=dudx; dvdy=dudx;
% end

t4 = 1/(nu^2-1);
t6 = 1/(1+nu);
t7 = E*dudy.*t4*.25;
t8 = dudx+1.0;
t9 = E*dudy.*t4.*t8*.25;
t10 = dvdy+1.0;
t11 = t7+t9-E*dvdx.*t6.*t10*.25;
t12 = dvdy.*2.0;
t13 = dudx.^2;
t14 = dudy.^2;
t15 = dvdy.^2;
t16 = dvdx.^2;
t17 = E*dvdx.*t4.*(1.0./2.0);
t18 = E*dudx.*dvdx.*t4*.25;
t19 = t17+t18-E*dudy.*t6*.25-E*dudy.*dvdy.*t6.*(1.0./8.0);
t20 = E*dudy.*dvdx.*nu.*t4*.25;
t21 = t20-E*t6.*(1.0./2.0)-E*dudx.*t6*.25-E*dvdy.*t6*.25-E*dudx.*dvdy.*t6.*(1.0./8.0);
t22 = dudx.*2.0;
t23 = dvdy+2.0;
t24 = nu-1.0;
t25 = E*nu.*t4;
t26 = E*dudx.*nu.*t4.*(1.0./2.0);
t27 = E*dvdy.*nu.*t4.*(1.0./2.0);
t28 = E*dudx.*dvdy.*nu.*t4*.25;
t29 = t25+t26+t27+t28-E*dudy.*dvdx.*t6.*(1.0./8.0);
t30 = E*dudy.*t4.*t23.*(1.0./8.0);
t31 = E*dudy.*dvdy.*t4.*(1.0./8.0);
t32 = t7+t30+t31-E*dvdx.*t6.*(1.0./8.0)-E*dvdx.*t4.*t8.*t24.*(1.0./8.0);
t33 = dudy.*2.0;
t34 = dvdx.*2.0;
t35 = dudx.*dudy.*2.0;
t36 = dvdx.*dvdy;
t37 = t33+t34+t35+t36;
t38 = 1.0./t24;
t39 = E*dvdx.*t23.*t38.*(1.0./8.0);
t40 = t39-E*t6.*t37.*(1.0./8.0);
out1 = [E*t4.*(dudx.*6.0+t13.*2.0+t14+t16+4.0)*.25+E*nu.*t4.*(t12+t15)*.25;
  t11;
  t19;
  t29;
  t11;
  E*t4.*(t12+t13+t14.*2.0+t15+t22+2.0)*.25+E*nu.*t4.*(t16-2.0)*.25;
  t21;
  t32;
  t19;
  t21;
E*t4.*(t12+t13+t15+t16.*2.0+t22+2.0)*.25+E*nu.*t4.*(t14-2.0)*.25;
t40;
t29;
t32;
t40;
E*t4.*(dvdy.*6.0+t14+t15.*2.0+t16+4.0)*.25+E*nu.*t4.*(t13+t22)*.25];

c = -out1([1 5 6 9 10 13 14 11 15 16], :);

end