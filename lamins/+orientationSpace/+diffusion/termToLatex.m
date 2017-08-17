function [ latex ] = termToLatex( s )
%termToLatex Summary of this function goes here
%   Detailed explanation goes here

import orientationSpace.diffusion.*;

if(~isscalar(s))
    latex = arrayfun(@termToLatex,s(:),'UniformOutput',false);
    return;
end

latex = '';

if(s.coeff ~= 1)
    latex = [latex num2str(s.coeff) ' '];
end

if(s.D == 1)
    latex = [latex 'D '];
elseif(s.D > 1)
    latex = [latex 'D^' num2str(s.D) ' '];
end

for i=1:length(s.rhod)
    if(i == 1)
        if(s.rhod(i) == 1)
            latex = [latex sprintf('\\pdv{\\rho}{\\theta} ')];
        elseif(s.rhod(i) > 1)
            latex = [latex sprintf('\\left(\\pdv{\\rho}{\\theta}\\right)^%d ',s.rhod(i))];
        end
    else
        if(s.rhod(i) == 1)
            latex = [latex sprintf('\\pdv[%d]{\\rho}{\\theta} ',i)];
        elseif(s.rhod(i) > 1)
            latex = [latex sprintf('\\left(\\pdv[%d]{\\rho}{\\theta}\\right)^%d ',i,s.rhod(i))];
        end
    end
end

for i=1:length(s.timed)
    if(i == 1)
        if(s.timed(i) == 1)
            latex = [latex sprintf('\\dv{t}{\\theta} ')];
        elseif(s.timed(i) > 1)
            latex = [latex sprintf('\\left(\\dv{t}{\\theta}\\right)^%d ',s.timed(i))];
        end
    else
        if(s.timed(i) == 1)
            latex = [latex sprintf('\\dv[%d]{t}{\\theta} ',i)];
        elseif(s.timed(i) > 1)
            latex = [latex sprintf('\\left(\\dv[%d]{t}{\\theta}\\right)^%d ',i,s.timed(i))];
        end
    end
end

end

