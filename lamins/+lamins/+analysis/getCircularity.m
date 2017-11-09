function [ c ] = getCircularity( MD )
%getCircularity Obtains the circularity score for all the data in
%ML

if(~isa(MD,'MovieData'))
    ML = lamins.util.toMovieList(MD);
    c = cellfun(@lamins.analysis.getCircularity,ML.movies_,'Unif',false);
    c = [c{:}];
    return;
end

skeletons = load([MD.outputDirectory_ filesep 'skeletons_2015_06_10.mat']);

% expand tz locations by channel
ctz = lamins.util.TZtoCTZ_LinearInd(MD,skeletons.tz);

% get circularity from mask
L = lamins.classes.LaminsData(MD);
images = L.getImages;
images = images(ctz);
c.maskCircularity = zeros(size(images));
c.maskCircularity = arrayfun(@(x) x.getMaskCircularity,images,'ErrorHandler',@(varargin) NaN);

% get circularity from final skeletons
S = skeletons.S3(ctz);
c.skelCircularity = zeros(size(S));
c.skelCircularity = cellfun(@(x) x.getNuclearCircularity,S,'ErrorHandler',@(varargin) NaN);

end

