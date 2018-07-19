% Square matrix, where Cost(i,j) is the cost of classifying a point into class j if its true class is i 
% (i.e., the rows correspond to the true class and the columns correspond to the predicted class). 
% To specify the class order for the corresponding rows and columns of Cost, additionally specify 
% the ClassNames name-value pair argument.
function costVal = LCH_getCostVal(n1,n2)
n = n1 + n2;
costVal = [...
    [0,         1-n1/n];
    [1-n2/n,    0];
    ];
end
