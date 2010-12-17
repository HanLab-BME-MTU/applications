function [center,fx_sum,fy_sum,ftot_sum,sx_sum,sy_sum,nVec_mean,fx_int,fy_int,ftot_int,nVec]=calcIntfacialStress(curve,sxx,syy,sxy,p,pixSize_mu,method)
% input:    curve has to be a nx2 matrix, first column beeing the x-component,
%           second column being the y-component.
% center:   center of each line segment, where the interfacial force is
%           located.

numPieces=size(curve,1)-1;

% The center of each line segment
center=1/2*(curve(2:end,:)+curve(1:end-1,:));

% The connection vector for each line segment:
segVecs=curve(2:end,:)-curve(1:end-1,:);

% calculate the length of each line segments:
dl=sqrt(sum(segVecs.^2,2));

% calculate normal vectors for each line segment:
x1=-segVecs(:,2);
x2= segVecs(:,1);
nFac=sqrt((x1.^2+x2.^2));

% Now go and normalize the vectors:
x1=x1./nFac;
x2=x2./nFac;
nVec=horzcat(x1,x2);
nVec_mean=mean(nVec);

% figure(10)
% plot(curve(:,1),curve(:,2),'-k')
% hold on
% quiver(center(:,1),center(:,2),nVec(:,1),nVec(:,2))
% hold off

% interpolate the stresses at each center of the line segments:

sxx_func = TriScatteredInterp(p(1,:)',p(2,:)',sxx,method);
syy_func = TriScatteredInterp(p(1,:)',p(2,:)',syy,method);
sxy_func = TriScatteredInterp(p(1,:)',p(2,:)',sxy,method);

isxx=sxx_func(center(:,1),center(:,2));
isyy=syy_func(center(:,1),center(:,2));
isxy=sxy_func(center(:,1),center(:,2));


% Calculate the x- and y-component of the force on each line segement using 
% a very simple discretization of the integral.
% for segment k it would be: 
% fx(k)=sxx(k)*n(k)*l(k)+sxy(k)*n(k)*l(k)
% fx(k)=sxx(k)*n(k)*l(k)+sxy(k)*n(k)*l(k)
% In vector notation:
fx_sum=-(isxx.*nVec(:,1).*dl+isxy.*nVec(:,2).*dl);
fy_sum=-(isxy.*nVec(:,1).*dl+isyy.*nVec(:,2).*dl); % using syx=sxy
ftot_sum=[sum(fx_sum) sum(fy_sum)];
% the stress (force normalized by the segment length) is:
sx_sum=fx_sum./dl;
sy_sum=fy_sum./dl;

% !!! The following has to be calculated after calculating sx_sum and sy_sum
% in the line above!: 
% fx_sum, fy_sum, ftot_sum have the dimension of a force, but still lack
% the factor: 
factor=pixSize_mu^2*10^(-3);
% to bring them to the unit nN. That is done now:
fx_sum=fx_sum*factor;
fy_sum=fy_sum*factor;
ftot_sum=ftot_sum*factor;

% Instead of the linearization above we could also use quad to explicitely
% calculate the integral:
fx_int=[];
fy_int=[];
ftot_int=[];
