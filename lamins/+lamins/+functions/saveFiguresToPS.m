function [ out ] = saveFiguresToPS( ML )
%saveFiguresToPS Saves multiple figures to a single PS file
%
% ML is a MovieList, or string for MovieList location

if(ischar(ML))
    ML = MovieList.load(ML);
end

if(isa(ML,'MovieList'))

    if(isempty(ML.movies_))
        ML.sanityCheck();
    end
    
end

if(isa(ML,'MovieData'))
    % If MovieData, then created a fake struct and put the Movie in a cell
    % array
    ML = struct('movies_',{num2cell(ML)},'movieListPath_',ML.outputDirectory_);
end

analysisDate = '2015_06_10';

import lamins.functions.*;

tzs = getTZ(ML);
pathComp = strsplit(ML.movieListPath_,filesep);
file = [ML.movieListPath_ filesep pathComp{end} '.ps'];
pdffile = [ML.movieListPath_ filesep pathComp{end} '.pdf'];
delete(file);

for i=1:length(ML.movies_)
    MD = ML.movies_{i};
    for ch = 1:length(MD.channels_)
        linidx = sub2ind([length(MD.channels_) MD.nFrames_ MD.zSize_],ch,tzs(i));
        h = openfig([MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(linidx) '.fig']);
        set(h,'Units','inches');
        set(h,'Position',[0 0 10.5 8]);
        set(h,'PaperOrientation','landscape');
        set(h,'PaperPositionMode','manual');
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[0.25 0.25 10.5 8]);
        text(50,110,['CTZ Frame/Slice ' num2str(linidx)],'Color','y','Interpreter','none');
        drawnow;
        print(h,file,'-dpsc','-append');
        close(h);
    end
end

system(['ps2pdf ' file ' ' pdffile])

out = tzs;

end

