function [filtered]=filterMany(Img,Masks)
%filterMany, applies multiple filters to the same image (Img) and adds the
%responses, filtered.
%
%Inputs : Img -> a 2D matrix
%         Masks -> a cell array containing multiple masks
%
%Output: filtered -> sum of all filtering responses
%
%
%Jeffrey Werbin
%Harvard Medical Schoool
%
% last update 12/13/12

n = numel(Masks);

filtered = zeros(size(Img));

for i= 1:n
    temp = filter2(Masks{i},Img);
    filtered = max(filtered,temp);
end


end