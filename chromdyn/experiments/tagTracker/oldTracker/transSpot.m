function  [model]=transSpot(match,dataProperties)


% Constants
%global FT_SIGMA;
if nargin <2 | isempty(dataProperties)
    %load_sp_constants;
    error('Missing dataProperties');
else
    FT_SIGMA = dataProperties.FT_SIGMA;
end


model=match;


% compute new coords


for s=1:length(model)
     model(s).center=match(s).center+match(s).parms;
     model(s).coords=match(s).coords+ones(length(match(s).coords),1)*match(s).parms;
end;


% forward map of gaussians with current parms to calc ratio in template

%  calc for overlaps
for s=1:length(match)
    %tempVal=match(s).gauss;
    model(s).gauss=gauss_3D(model(s).amp,FT_SIGMA,model(s).center,model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3));
    tempVal=model(s).gauss;
    for ov=setdiff(1:length(match),s)
        gval=gauss_3D(model(ov).amp,FT_SIGMA,model(ov).center,model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3));
        tempVal=tempVal+gval;
    end;
    model(s).ratio=model(s).gauss./tempVal;
end;



% compute gradient

