function [ code ] = termsToCode( s, partialName, fullName )
%termsToCode Convert derivative terms to code

if(nargin < 2)
    partialName = 'rhodv';
end
if(nargin < 3)
    if(isfield(s,'timed'))
        fullName = 'dnt_dmn';
    else
        fullName = 'dnm_dtn';
    end
end


import orientationSpace.diffusion.*;

if(~isscalar(s))
    code = arrayfun(@(s) termsToCode(s,partialName,fullName),s(:),'UniformOutput',false);
    if(nargout == 0)
        disp(['  ' strjoin(code,' ...\n+ ')]);
    end
    return;
end

code = '';

if(s.coeff ~= 1)
    code = [code sprintf('%2d',s.coeff) ' .* '];
else
    code = [code '      '];
end

if(s.D == 1)
    code = [code 'D   .* '];
elseif(s.D > 1)
    code = [code 'D^' num2str(s.D) ' .* '];
else
    code = [code '       '];
end

for i=1:length(s.rhod)
    if(s.rhod(i) == 1)
        code = [code sprintf('%s(:,:,%d) .* ',partialName,i)];
    elseif(s.rhod(i) > 1)
        code = [code sprintf('%s(:,:,%d).^%d .* ',partialName,i,s.rhod(i))];
    end
end

if(~isfield(s,'timed'))
    s.timed = [];
end

suffix = ' .* ';
for i=1:length(s.timed)
%     if(i == length(s.timed))
%         suffix = '';
%     end
    if(s.timed(i) == 1)
        code = [code sprintf('%s(:,:,%d)   %s',fullName,i,suffix)];
    elseif(s.timed(i) > 1)
        code = [code sprintf('%s(:,:,%d).^%d%s',fullName,i,s.timed(i),suffix)];
    end
end
if(~isempty(s.timed))
    code = code(1:end-length(suffix));
end

if(~isfield(s,'thetad'))
    s.thetad = [];
end

suffix = ' .* ';
for i=1:length(s.thetad)
%     if(i == length(s.thetad))
%         suffix = '';
%     end
    if(s.thetad(i) == 1)
        code = [code sprintf('%s(:,:,%d)   %s',fullName,i,suffix)];
    elseif(s.thetad(i) > 1)
        code = [code sprintf('%s(:,:,%d).^%d%s',fullName,i,s.thetad(i),suffix)];
    end
end
if(~isempty(s.thetad))
    code = code(1:end-length(suffix));
end

    if(nargout == 0)
        disp(code);
    end

end

