% Square matrix, where Cost(i,j) is the cost of classifying a point into class j if its true class is i 
% (i.e., the rows correspond to the true class and the columns correspond to the predicted class). 
% To specify the class order for the corresponding rows and columns of Cost, additionally specify 
% the ClassNames name-value pair argument.
function costVal = LCH_getCostValMultLabels(ns)
n = sum(ns);
nlabels = length(ns);

costVal = nan(nlabels);

for i = 1 : nlabels
    costVal(i,:) = ones(1,nlabels).* (1-ns(i)/n);
    costVal(i,i) = 0;
end
end
