load('movieLists.mat');
import lamins.functions.*;
tables = cell(1,length(MLs));
for i=1:32
    tables{i}.tz = getTZ(MLs(i));
    try
        tables{i}.tz2 = getTZ(MLs(i),2);
    catch err
        tables{i}.tz2 = Inf(size(tables{i},1),1);
    end
    tables{i}.mintz = min(tables{i}.tz,tables{i}.tz2);
end

T = vertcat(tables{:});
save('tzTable.mat','T','tables');
writetable(T,'tzTable.csv','QuoteStrings',true);