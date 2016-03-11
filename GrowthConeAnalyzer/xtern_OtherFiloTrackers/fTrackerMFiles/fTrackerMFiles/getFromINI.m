function outWithFormat=getFromINI(readCell, varName)

% Reads variables from the cell obtained after using INIfile. Lookes for
% the variable name and assings teh avlue according to the desired type.
% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


indice=strmatch(lower(varName), readCell, 'exact')-size(readCell,1)*2;

if size(str2num(readCell{indice,4}),1)==0
    outWithFormat=readCell{indice,4};
else
    outWithFormat=str2num(readCell{indice,4});
end

% switch varType
%     case 'text'
%         outWithFormat=readCell{indice,4};
%     case 'number'
%         outWithFormat=str2num(readCell{indice,4});
% end
       
    

