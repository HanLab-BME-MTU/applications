function [data,labels]=makeFeaturesHighOrder(data,labels,combos)
% takes some data (n x m), labels (1 x m) and combos ( k x m )
% uses combos to indicate how the columns of data and labels should be combined to make new features
% 


[n,m] = size(data);
k = size(combos,1);

for i=1:k
	temp = combos(i,:);
        [~,j]=max(temp);
        nFeat = data(:,j);
        nLabel = labels{j};        
        temp(j) = temp(j)-1;
        while(sum(temp)>0)
          [~,j]=max(temp);
          nFeat = nFeat.*data(:,j);
          nLabel = [nLabel,'*',labels{j}];
          temp(j) = temp(j)-1;
        end
	data = [data,nFeat];
        labels = [labels,{nLabel}];
end

