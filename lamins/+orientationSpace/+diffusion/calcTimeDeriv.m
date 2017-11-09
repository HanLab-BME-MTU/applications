function [ s ] = calcTimeDeriv( q, n )
%calcTimeDeriv Calculate equation for full derivative with respect to
%time of order q of the partial derivative with respect to theta of order
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
%  .thetad - theta derivative portion, number identifies the exponent, column
%           is the order of the derivative of time with respect to theta
%  .D     - number describes the exponent of the diffusion coefficient
%
% Example OUTPUT
%
%% First full derivative of the first partial derivative
% s = calcTimeDeriv(1,1);
%% \pdv[2]{\rho}{\theta} \dv{\theta}{t}
% s(1)
% 
% ans = 
% 
%       rhod: [0 1]
%      coeff: 1
%     thetad: 1
%          D: 0
% 
%% D \pdv[3]{\rho}{\theta} 
% s(2)
% 
% ans = 
% 
%       rhod: [0 0 1]
%      coeff: 1
%     thetad: 0
%          D: 1
%% Overall \pdv[2]{\rho}{\theta} \dv{\theta}{t} + D \pdv[3]{\rho}{\theta} 
% 


import orientationSpace.diffusion.*;

if(q == 0)
    % 0th order derivative is just the nth partial derivative
    s.rhod(1,n) = 1;
    s.coeff(1,1) = 1;
    s.thetad(1,1) = 0;
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
    s = calcTimeDeriv(q-1,n);
    ss = cell(1,length(s));
    % differentiate each term
    for i=1:length(s)
        ss{i} = sdiff(s(i));
    end
    s = [ss{:}];
end
    
%% Simplify and sort
s = normalize(s);
s = simplify(s);
try
    s = sortTerms(s);
catch err
    % sorting is auxillary
end

end

function s = sdiff(s)
    % find non-zero t derivatives with respect to theta
    nztd = s.thetad ~= 0;
    if(~any(nztd)) % no t derivatives with respect to theta
        % partial deriv with respect to theta
        s(1).rhod = [0 s(1).rhod];
        % partial deriv with respect to time
        s(2) = s(1);
        s(2).rhod = [0 s(2).rhod];
        s(2).D = s(2).D + 1;
        % add derivative of theta with respect to time
        s(1).thetad = 1;
    else
        ss = s(1);
        % differentiate the rho derivative portion first
        thetad = ss.thetad;
        ss.thetad = 0;
        ss = sdiff(ss);
        % multiply by the theta derivative portion
        for i=1:length(ss)
            l = max(length(ss(i).thetad),length(thetad));
            ss(i).thetad(end+1:l) = 0;
            thetad(end+1:l) = 0;
            ss(i).thetad = ss(i).thetad+thetad;
        end
        % differentiate the theta derivative portion
        td_idx = find(nztd);
        sss(length(td_idx)) = s(1);
        s(1).thetad(end+1) = 0;
        template = zeros(1,length(s(1).thetad));
        for i=1:length(td_idx)
            sss(i) = s(1);
            ex = s(1).thetad(td_idx(i));
            sss(i).thetad = template;
            sss(i).coeff = sss(i).coeff * ex;
            sss(i).thetad(td_idx(i)) = -1;
            sss(i).thetad(td_idx(i)+1) = 1;
            sss(i).thetad = sss(i).thetad + s(1).thetad;
        end
        s = [ss sss];
    end
end

function s = normalize(s)
    max_td_length = max(arrayfun(@(s) length(s.thetad),s));
    for i=1:length(s)
        s(i).thetad(end+1:max_td_length) = 0;
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

