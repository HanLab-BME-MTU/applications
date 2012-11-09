function [nx_n, ny_n]=prCalculateNormal(img, sp_x_n, sp_y_n, r_n)
% prCalculateNormal calculates the unit normals relative to the edge 
%               This function calculates the unit normals along a edge
%               given by the two component spline at the paremeters given
%               by "r_n". The normal points away from the object assuming
%               that the object has a higher average intensity.
%
% SYNOPSIS    [nx_n, ny_n, band_mask]=prGetProtRegion(img, sp_x_n, sp_y_n, r_n)
%
% INPUT       img       : normalized grayscale image
%             sp_x_n    : edge x B-spline 
%             sp_x_n    : edge y B-spline 
%             r_n       : spline parameter
%
% OUTPUT      nx_n      : x-component of the edge unit normal
%             ny_n      : y-component of the edge unit normal
%   
% First Created by Matthias Machacek 11/11/03 (prGetProtRegion.m)
% Modefied as prCalculateNormal.m by Shann-Ching Chen, 05/22/2008

[n_img, m_img]=size(img);
knots_nr=length(r_n);

x_n=fnval(sp_x_n,r_n);
y_n=fnval(sp_y_n,r_n);

%derivative spline 
sp_dx_n = fnder(sp_x_n);
sp_dy_n = fnder(sp_y_n);
%derivatives at discrite locations
dx_n=fnval(sp_dx_n,r_n);
dy_n=fnval(sp_dy_n,r_n);
%normalize
l=sqrt(dx_n.^2+dy_n.^2);
dx_nn=dx_n./l;
dy_nn=dy_n./l;
%the normal unit vector (not oriented!!)
nx_n= dy_n;
ny_n=-dx_n;

%determine the object side of the edge
p1_x=round(x_n+2*nx_n);
p1_y=round(y_n+2*ny_n);

p2_x=round(x_n-2*nx_n);
p2_y=round(y_n-2*ny_n);

in_out1=0;
in_out2=0;

for i=1:length(p1_x)
    logic_exp1 = p1_x(i) >= 1     & p1_y(i) >= 1     & p2_x(i) >= 1     & p2_y(i) >= 1;
    logic_exp2 = p1_x(i) <= m_img & p1_y(i) <= n_img & p2_x(i) <= m_img & p2_y(i) <= n_img;
    if logic_exp1 & logic_exp2
        in_out1 = in_out1 + img(p1_y(i), p1_x(i));
        in_out2 = in_out2 + img(p2_y(i), p2_x(i));
    end
end

%assume that the background has a lower intensity
%the normal is pointin away from the object!
if in_out1 > in_out2
   nx_n= -nx_n;
   ny_n= -ny_n; 
end
