function [metaData] = pcMetaData()
metaData.dates = {'131121','131122','131123','131203','131204','131223','131225'};

metaData.day1 = {'131121','131203'};
metaData.day2 = {'131122','131204'};
metaData.day3 = {'131123'};

metaData.tumors = {'m481','m214','m528','m610','um12', 'm514', 'm405'}; % m528 is m530

metaData.highMetPot = {'m481','m214','um12', 'm514', 'm405'};
metaData.lowMetPot = {'m528','m610'};

end