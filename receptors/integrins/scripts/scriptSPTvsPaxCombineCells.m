
for i = 1 : length(movieStructAlphaVPax);
    
    disp(num2str(i))
    
    if movieStructAlphaVPax(i).activityLevel > 0
        
        tmp = movieStructAlphaVPax(i).fileName{1};
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisFAs'])
        
        cellRes(i,1) = load('sptDiffVsPaxIntRes.mat');

    end
    
end
