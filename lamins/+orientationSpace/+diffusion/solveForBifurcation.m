function [ bifurcation_orientation, bifurcation_K, bsol,fval,exitflag,output,jacobian ] = solveForBifurcation( response, Korg, lm, K, freq )
%solveForBifurcation Solves for the bifurcation point 

if(nargin < 3)
    lm = interpft_extrema(response);
    K = Korg;
end
if(nargin < 4)
    K = Korg;
end
if(nargin < 5)
    freq = false;
end

if(isempty(K))
    K = Korg;
end

K = repmat(K,size(lm)./size(K));

validLM = ~isnan(lm);
t = 1./(2*K(validLM)+1).^2;

init = [lm(validLM).'; t(:).'];
opt = optimoptions('fsolve');
opt.SpecifyObjectiveGradient = true;

f = orientationSpace.diffusion.generateBifurcationFunction(response,Korg,freq);

[bsol,fval,exitflag,output,jacobian] = fsolve(f,init,opt);

bsol(:,any(abs(fval) > opt.FunctionTolerance)) = NaN;

bifurcation_orientation = lm;
bifurcation_orientation(validLM) = bsol(1,:);
bifurcation_K = lm;
bifurcation_K(validLM) = (1./sqrt(bsol(2,:))-1)/2;

end

