function [wt]=multiWeibullODF(t,params)

% multiple Weibull original distribution function

offset = params(1);
% what's the role of the offset? if the offset is the offset in the
% cumulative distribution function, then it translates into differences in
% the normalization in the original curve??
% Because of this issue, set offset to zero before the function is carries
% out if necessary

n = round((length(params)-1)/3);

for w=1:n
    % the original probability distribution function for a
    % general Weibull is, analytically, 
    % (k/l)*(x/l)^(k-1)*exp(-x/lambda)^k
    % This function is normalized, i.e. the area under the function
    % integrated from 0 to inf equals 1 - thus, the amplitude amp is the
    % area scaling factor
    
    amp(w) = params(1+3*w-2);
    lambda(w) = params(1+3*w-1);
    ka(w) = params(1+3*w);
                        
    weibullMat(w,:) = amp(w) * (ka(w)/lambda(w)) * ((t/lambda(w)).^(ka(w)-1)) .*...
        exp(-(t/lambda(w)).^ka(w));
end

wt = offset + sum(weibullMat,1);


end % of function