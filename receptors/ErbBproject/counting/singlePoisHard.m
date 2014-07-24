function y = singlePoisHard(t,indic,Pois)
% given a matrix Pois that contains the likelihood of the next event as a
% funciton of t-ti where the units are a standardized time step
% 
% indic is list of data points assigned to this model (varying between 0
%           and 1
%
% the output y(i) is Pois(t(i)-t(j)) where tj is the most recent apperance
% before time ti
%
tmod = t(indic);
[dim, nmod] = size(tmod);
[dim, npoints] = size(t);

[h, bin] = hist(indic.*t,[1:max(t)]); % this gives the number of points a time t in this model
ch = cumsum(h);

if dim ~= 1
    error('time should only have one dimension');
end

%computs all delta t between all point and points in model
dt = repmat(t,[nmod,1])-repmat(tmod',[1,npoints]);
dt(dt<0)=NaN;
%this ensures that if the point is assigned to the model that it is
%evaluated based on the nearest point that came before it.

%This avoids a point assigned to the model being considered a second point
%at time t.
dt(dt==0 & repmat(indic,[nmod,1]))=NaN; 

%Finds the smallest value of dt >= 0 
%if best == zero means another point appeared at this time and they can't
%belong to the same cluster
best = min(dt);
best(isnan(best))=0;

%calulates probability, adds 1 to handle the 0 index
y = Pois(best+1);


%special case first point, penalizes being the first point in the cluster
%your probability will always be zero
% To fix this 
%if a point i that has a y of 0 and their are no points in the model that 
%came before point i then y is set to 1 

y(y==0 & (ch(t)==0))=1;


end



