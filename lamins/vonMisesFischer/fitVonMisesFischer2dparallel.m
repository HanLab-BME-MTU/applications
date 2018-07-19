function [ fits ] = fitVonMisesFischer2dparallel( data, alpha, maxN )
%fitVonMisesFischer2dparallel Do fits in parallel

% each row is fit and a cell array of outputs is returned

data = num2cell(data,2);
data = distributed(data);
fits = cellfun(@(d) fitVonMisesFischer2d(d,alpha,maxN),data,'UniformOutput',false);

fits = gather(fits);


end

