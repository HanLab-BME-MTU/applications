function[data,labels] = makeExtraFeatures(data,labels)
% data is a m x n matrix (with n features)
% labels is 1 x n cell array of string with the names of the feautres
%
% returns both back with all 2nd order features appened
%

  n = size(data,2);
n1 = size(labels,2);

if n ~= n1
  error('data and labels need to have the same number of columns');
end

for i = 1:n
	  for j=i:n
		  data = [data,data(:,i).*data(:,j)];
                  labels = [labels,{[labels{i},'*',labels{j}]}];
          end
end

end
