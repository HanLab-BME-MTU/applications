function [ s ] = calcThetaDeriv( q, n )
%calcThetaDeriv Calculate equation for full derivative with respect to
%theta of order q of the partial derivative with respect to theta of order
%n for the heat/diffusion equation
%
% INPUT
% q - order of the full derivative
% n - order of the partial derivative
%
% OUTPUT
% s - structure array describing each additive term
%  .rhod - number identifies the exponent, column is order of partial
%          derivative of the response with respect to theta
%  .coeff - scalar multiplicative factor
%  .timed - time derivative portion, number identifies the exponent, column
%           is the order of the derivative of time with respect to theta
%  .D     - number describes the exponent of the diffusion coefficient
%
% Example OUTPUT
%
%% First full derivative of the first partial derivative
% s = calcThetaDeriv(1,1);
%% \pdv[2]{\rho}{\theta}
% s(1)
% 
% ans = 
% 
%      rhod: [0 1]
%     coeff: 1
%     timed: 0
%         D: 0
% 
%% D \pdv[3]{\rho}{\theta} \dv{t}{\theta}
% s(2)
% 
% ans = 
% 
%      rhod: [0 0 1]
%     coeff: 1
%     timed: 1
%         D: 1
%% Overall \pdv[2]{\rho}{\theta} + D \pdv[3]{\rho}{\theta} \dv{t}{\theta}
% 


import orientationSpace.diffusion.*;

if(q == 0)
    % 0th order derivative is just the nth partial derivative
    s.rhod(1,n) = 1;
    s.coeff(1,1) = 1;
    s.timed(1,1) = 0;
    s.D(1,1) = 0;
    return;
% elseif(q == 1)
%     s(1) = calcThetaDeriv(q-1,n+1);
%     
%     s(2) = calcThetaDeriv(q-1,n+2);
%     s(2).D = s(2).D + 1;
%     s(2).timed(:,1) = s(2).timed(:,1) + 1;
else
    % Calculate the lower order derivative
    s = calcThetaDeriv(q-1,n);
    ss = cell(1,length(s));
    % differentiate each term
    for i=1:length(s)
        ss{i} = sdiff(s(i));
    end
    s = [ss{:}];
end
    
%% Simplify and sort
s = simplify(s);
try
    s = sortTerms(s);
catch err
    % sorting is auxillary
end

end

function s = sdiff(s)
    % find non-zero t derivatives with respect to theta
    nztd = s.timed ~= 0;
    if(isscalar(s.timed) || isscalar(s.timed(nztd)))
        if(all(s.timed == 0)) % no t derivatives with respect to theta
            % partial deriv with respect to theta
            s(1).rhod = [0 s(1).rhod];
            % partial deriv with respect to time
            s(2) = s(1);
            s(2).rhod = [0 s(2).rhod];
            s(2).D = s(2).D + 1;
            s(2).timed(:,1) = s(2).timed(:,1)+1;
        else % one t derivative with respect to theta
            ss = s(1);
            % differentiate the rho derivative portion first
            timed = ss.timed;
            ss.timed = 0;
            ss = sdiff(ss);
            % multiply by the t derivative portion
            for i=1:length(ss)
                ss(i).timed = ss(i).timed+timed;
            end
            % differentiate the t derivative portion
            sss = s(1);
            if(sss.timed(nztd) == 1)
                % no exponent, just bump up the derivative
                sss.timed = [0 sss.timed];
            else
                % exponent
                ex = sss.timed(nztd);
                % multiply by exponent
                sss.coeff = sss.coeff * ex;
                % decrease exponent by one
                sss.timed(nztd) = sss.timed(nztd)-1;
                % multiply by next derivative
                sss.timed([false nztd]) = 1;
            end
            s = [ss sss];
        end
    else
        % TODO: compound time derivative
        keyboard;
    end
end

function ss =  simplify(ss)
    % Combine similar additive terms
    if(length(ss) < 2)
        return;
    end
    % Remove scalar coefficient
    temp = rmfield(ss,'coeff');
    same = arrayfun(@(s) isequal(temp(1),s),temp(2:end));
    % Combined term now has the sum of the scalar coefficients
    ss(1).coeff = sum([ss([true same]).coeff]);
    % Remove extra terms
    ss([false same]) = [];
    % Recursive
    ss = [ss(1) simplify(ss(2:end))];
end

function s = sortTerms(s)
% only works if all the terms have only one partial derivative
    [~,idx] = sort(arrayfun(@(s) find(s.rhod),s));
    s = s(idx);
end

