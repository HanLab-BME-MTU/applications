function [prmVect, prmStd, y] = fitExp(t, signal, prmVect, prmSel, mode)

if nargin<5
    mode = '+';
end

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

estIdx = false(1,3); % [k A dt]
estIdx(regexp('kAt', ['[' prmSel ']'])) = true;

[x,resnorm,~,~,~,~,J] = lsqnonlin(@costExp, prmVect(estIdx), [], [], opts, signal, t, prmVect, estIdx, mode);
prmVect(estIdx) = x;
prmVect(1:3) = deal(abs(prmVect(1:3)));

C = resnorm*full(inv(J'*J));
prmStd = sqrt(diag(C)/(numel(signal)-length(x) - 1));

k = prmVect(1);
A = prmVect(2);
dt = prmVect(3);
if strcmp(mode, '+')
    y = A*(1 - exp(-k*(t-dt)));
else
    y = A*exp(-k*(t-dt));
end

function [v] = costExp(p, signal, t, prmVect, estIdx, mode)
prmVect(estIdx) = p;
k = abs(prmVect(1));
A = abs(prmVect(2));
dt = prmVect(3);
if strcmp(mode, '+')
    v = signal - A*(1-exp(-k*(t-dt)));
else
    v = signal - A*exp(-k*(t-dt));
end