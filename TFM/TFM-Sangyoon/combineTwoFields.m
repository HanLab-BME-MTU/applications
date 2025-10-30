function [displField] = combineTwoFields(displField1, displField2)
%function [displField] = combineTwoFields(displField1, displField2)
%combines two fields
% Sangyoon Han, June, 2020
numFields1 = numel(displField1);
numFields2 = numel(displField2);

if numFields1 ~= numFields2
    error('The two fields are different in number!')
else
    displField = displField1;
    for ii=1:numFields1
        displField(ii).pos = [displField1(ii).pos; displField2(ii).pos];
        displField(ii).vec = [displField1(ii).vec; displField2(ii).vec];
    end
end



