tsFiles={imdsValid.Files};
trFiles={imdsTrain.Files};
setTR=cellfun(@(x) x(end), cellfun(@(x) strsplit(x,'/'),trFiles{:},'Unif',0));
tsFiles={imdsValid.Files};
setTS=cellfun(@(x) x(end), cellfun(@(x) strsplit(x,'/'),tsFiles{:},'Unif',0));
intersect(setTR,setTS);