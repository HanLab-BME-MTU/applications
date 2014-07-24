%function [I,Istd]=corrFunc(u_rp,u_rp_r,meanSub)
function [I,Isem95,Istd]=corrFunc(u_rp_x,u_rp_r_x,u_rp_y,u_rp_r_y,meanSub)
% correlation function for one (mean-substracted) component of the velocity
% I(r)=<u(rp+r).u(rp+r)>/sqrt(<u^2(rp+r)> <u^2(rp+r)>);
doNorm=1;

if nargin <3 || (isempty(u_rp_y) &&  isempty(u_rp_r_y))
    if nargin > 4 && meanSub==1
        % calculate the mean substracted velocities:
        u_rp_x  =u_rp_x  -nanmean(u_rp_x);
        u_rp_r_x=u_rp_r_x-nanmean(u_rp_r_x);
    end
    
    if doNorm==0
        % calculate the correlation:
        I    = nanmean(u_rp_x.*u_rp_r_x);%/sqrt(nanmean(u_rp.^2) * nanmean(u_rp_r.^2));
        Istd =  nanstd(u_rp_x.*u_rp_r_x);%/sqrt(nanmean(u_rp.^2) * nanmean(u_rp_r.^2));
    else
        I    = nanmean(u_rp_x.*u_rp_r_x)/sqrt(nanmean(u_rp_x.^2) * nanmean(u_rp_r_x.^2));
        Istd =  nanstd(u_rp_x.*u_rp_r_x)/sqrt(nanmean(u_rp_x.^2) * nanmean(u_rp_r_x.^2));
    end
else
    % correlation function for 2D (mean-subtracted velocity) vectors:
    % I(r)=<u(rp+r) x u(rp+r)>/sqrt(<u^2(rp+r)> <u^2(rp+r)>);
    % I am not sure if x above means the skalar product (or dot product) which
    % yields only one correlation value 'I' or if it should be the tensor
    % product, such that we get four values for I namely Ixx,Ixy,Iyx,Iyy.
    % Here we chose skalar product!

    if nargin > 4 && meanSub==1
        % calculate the mean substracted velocities. Here it is also not clear if
        % the mean of each component should be substracted or if the combined mean 
        % of the two components should be substracted:
        % Here we substract each mean separately!
        u_rp_x  =u_rp_x  -nanmean(u_rp_x);
        u_rp_y  =u_rp_y  -nanmean(u_rp_y);
        u_rp_r_x=u_rp_r_x-nanmean(u_rp_r_x);
        u_rp_r_y=u_rp_r_y-nanmean(u_rp_r_y);
    end

    % calculate the correlation:
    if doNorm==0
        I    = nanmean(u_rp_x.*u_rp_r_x + u_rp_y.*u_rp_r_y);%/sqrt(nanmean(u_rp_x.^2+u_rp_y.^2) * nanmean(u_rp_r_x.^2+u_rp_r_y.^2));
        Istd =  nanstd(u_rp_x.*u_rp_r_x + u_rp_y.*u_rp_r_y);%/sqrt(nanmean(u_rp_x.^2+u_rp_y.^2) * nanmean(u_rp_r_x.^2+u_rp_r_y.^2));
    else
        I    = nanmean(u_rp_x.*u_rp_r_x + u_rp_y.*u_rp_r_y)/sqrt(nanmean(u_rp_x.^2+u_rp_y.^2) * nanmean(u_rp_r_x.^2+u_rp_r_y.^2));
        Istd =  nanstd(u_rp_x.*u_rp_r_x + u_rp_y.*u_rp_r_y)/sqrt(nanmean(u_rp_x.^2+u_rp_y.^2) * nanmean(u_rp_r_x.^2+u_rp_r_y.^2));
    end
end
Isem95=facSEMtoSEM95*Istd/sqrt(sum(~isnan(u_rp_x)));