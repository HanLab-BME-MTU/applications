function russianForm=getRussianForm(connectivityMatrix)


% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca

[row, col]=find(connectivityMatrix~=0);

for itConnexions=1:size(row,1)
    russianForm(itConnexions, 1)=row(itConnexions);
    russianForm(itConnexions, 2)=col(itConnexions);
    russianForm(itConnexions, 3)=connectivityMatrix(row(itConnexions), col(itConnexions));
end

