function [wt]=multiWeibullCDF(t,params)

% multiple Weibull cumulative distribution function

offset = params(1);
n = round((length(params)-1)/3);

for w=1:n
    % the cumulative probability distribution function for a
    % general Weibull is, analytically, 1-exp(-x/lambda)^k
    % where lambda is the scaling factor, and k is the shape factor
    % Since the original Weibull function is already normalized (as is
    % obvious from the 1-... in the CDF), the amplitude reflects the
    % contribution of each individual population
    % NOTE: to prevent physcially meaningless populations, only positive
    % amplitudes and positive time constants are allowed
    
    amp(w) = abs(params(1+3*w-2));
    lambda(w) = abs(params(1+3*w-1));
    ka(w) = params(1+3*w);
                        
    weibullMat(w,:) = amp(w) * ( 1-exp( -(t/lambda(w)).^ka(w) ) );
end

wt = offset + sum(weibullMat,1);

end % of function