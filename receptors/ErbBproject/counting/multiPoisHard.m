function y = multiPoisHard(t,indic,fuzzy,Pois)
% Evaluates a mixture of Pois in time
% given a matrix Pois that contains the likelihood of the next event as a
% funciton of t-ti where the units are a standardized time step
% 
% indic are the points associated with this model
%
% fuzzy is the likelihood that point i is associated with this model
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
% Pois terms gives the likelihood that i is the next event of j
% repmat(indic) term mutliples each pois by the probabilities  that point j
% is from that model
y = Norm.*Pois(dt).*repmat(indic',[1,npoints]);
y = sum(y);%sums over one dimension
N = sum(indic); 

if N > 0
    y=y/N;
else
    y=zeros(size(y));
end

%special case first point, penalizes being the first point in the cluster
%your probability will always be zero, if a point that has a y of 0 and an
%indic of > 0.3 (meaning a strong association with a model) y is set to 1

y(y==0 & fuzzy>0.3)=1;

end



