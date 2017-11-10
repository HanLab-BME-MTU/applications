D = dir('MEF*');
D = D([D.isdir]);
for d = 1:length(D)
    cd(D(d).name);
    MD_dir = dir('*.nd2');
    MDs = cellfun(@MovieData.load,{MD_dir.name},'Unif',false);
    ML = MovieList(MDs,pwd);
    ML.movieListPath_ = pwd;
    ML.movieListFileName_ = [D(d).name '_list.mat'];
    ML.save();
    cd('..');
end