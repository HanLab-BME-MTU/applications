% colocalizationScatterSim(R2, nObj, c)
%
% Inputs: 
%         R2   : correlation coefficient
%         nObj : number of CCPs
%         c    : noise coefficient

% Francois Aguet, July 2010

function colocalizationScatterSim(R2, nObj, c)

N = round((-7-12*sqrt(R2) - 5*R2)/(R2-1));


capacity = unidrnd(N*ones(1,nObj));
class1 = binornd(capacity, 0.5);
class2 = capacity-class1;

if nargin==3
    %C = zeros(2,2);
    %while round(C(1,2)^2*100)/100 ~= round(R2*100)/100
    
    class1 = poissrnd(c*class1)/c;
    class2 = poissrnd(c*class2)/c;
    idx = class1==0 | class2==0;
    class1(idx) = [];
    class2(idx) = [];
    %C = corrcoef(class1,class2);
    %end
end

C = corrcoef(class1,class2);
fprintf('R2 = %f\n', C(1,2)^2);

% ma = 1.2*max(max(class1)/mean(class1), max(class2)/mean(class2));
ma = 3.5;
colocalizationScatterPlot(class1/mean(class1), class2/mean(class2), [0 ma 0 ma], 'red', 'green');

C = cov(class1, class2, 1);
fprintf('Theoretical: mean = %.4f, var = %.4f, cov = %.4f\n', (N+1)/4, (N+1)*(N+5)/48, (N+1)*(N-7)/48);
fprintf('Simulation : mean = %.4f, var = %.4f, cov = %.4f\n', mean(class1), C(1,1), C(1,2));

%print('-depsc2', '-r300', ['sim_R=' num2str(C(1,2)^2/(C(1,1)*C(2,2)), '%.2f') '.eps']);