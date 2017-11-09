function [Ilbp] = tlbp(Ilbp1,Ilbp2)
Ilbp = nan(size(Ilbp1));
c = 0;
for i = 0 : 9
    for j = 0 : 9
        Ilbp(Ilbp1 == i & Ilbp2 == j) = c;
        c = c + 1;
    end
end
end