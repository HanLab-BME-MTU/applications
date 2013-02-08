function y = multiPois(t,indic,Pois)
% Evaluates a mixture of Pois in time
% given a matrix Pois that contains the likelihood of the next event as a
% funciton of t-ti where the units are a standardized time step
% 
% indic is the likelihood that a point is associated with this model
%
% the output y(i) is sum over j (indic(j)*Pois(t(i)-t(j)))/#(indic(j) > 0) where Pois = 0 for ti < tj 
%
[dim, npoints] = size(t);

if dim ~= 1
    error('time should only have one dimension');
end

%computs delta t
dt = repmat(t,[npoints,1])-repmat(t',[1,npoints]);
Norm = dt > 0;
dt(dt <= 0)=1;

% Norm term removes any terms with t <=0
% Pois terms gives the likelihood that i is the next event of model j
% repmat(indic) term mutliples each pois by the probabilities  that point i
% is from that model
y = Norm.*Pois(dt).*repmat(indic',[1,npoints]);
y = sum(y); %sums over one dimension

%special case first point, penalizes being the first point in the cluster
%your probability will always be zero, if a point that has a y of 0 and an
%indic of > 0.3 (meaning a strong association with a model) y is set to 1

y(y==0 & indic>0.3)=1;

end



