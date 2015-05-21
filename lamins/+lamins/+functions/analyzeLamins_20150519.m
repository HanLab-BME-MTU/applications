function [ ] = analyzeLamins_20150519( MD )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    import lamins.functions.*;
    import lamins.classes.*;

    if(isa(MD,'MovieList'))
        ML = MD;
        for mdi = 1:length(ML.movies_)
            MD = ML.getMovie(mdi);
            analyzeLamins_20150519(MD);
        end
        return;
    else
        R = CellReader(MD.getReader());
        means = cellfun(@(x) mean(x(:)),R.toCell);
        [pks,loc] = findpeaks(means(:));
       
        L = LaminsData(MD);
        L.params.steerable.sigma = 5;
        images = L.getImages;
        S = cell(numel(images),1);
        S2 = cell(numel(images),1);
        for ii=loc'
            I = images(ii);
            S{ii} = I.skeleton;
            S{ii}.cleanup
            score = lamins.functions.scoreEdges(S{ii},I.flattenIntensity);
            S{ii}.deleteEdges(score < 0);
            
            S2{ii} = S{ii}.copy;
            thresh = I.maskThresh(double(I));
            S2{ii}.auditEdges(double(I),[],thresh,thresh/2);
            
            h = figure;
            imshow(I);
            S{ii}.drawEdgesAsLines([],'r');
            S2{ii}.drawEdgesAsLines([],'g');
            title(MD.movieDataFileName_);
            text(50,50,[MD.movieDataFileName_ ' analyzed on 2015_05_19'],'Color','y','Interpreter','none')
            savefig(h,[MD.outputDirectory_ filesep 'skeleton_20150519_' num2str(ii) '.fig']);
            saveas(h, [MD.outputDirectory_ filesep 'skeleton_20150519_' num2str(ii) '.pdf']);
            saveas(h, [MD.outputDirectory_ filesep 'skeleton_20150519_' num2str(ii) '.png']);
            close(h);
         
       end
       file = [MD.outputDirectory_ filesep 'skeletons_20150519.mat'];
       save(file,'S','S2');
    end

end

