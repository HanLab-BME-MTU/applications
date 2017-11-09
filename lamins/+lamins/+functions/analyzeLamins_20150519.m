function [ ] = analyzeLamins_20150519( MD )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    import lamins.functions.*;
    import lamins.classes.*;

    if(isa(MD,'MovieList'))
        ML = MD;
        for mdi = 1:length(ML.movies_)
            MD = ML.getMovie(mdi);
            try
                analyzeLamins_20150519(MD);
            catch err
                warning(['Failure on : ' MD.movieDataFileName_]);
                disp(err);
            end
        end
        return;
    else
        R = CellReader(MD.getReader());
        C = R.toCell;
        C = shiftdim(C,2);
        means = cellfun(@(x) mean(x(:)),C(:));
        [pks,loc] = findpeaks(means(:));
%         loc = 1;
        if(isempty(pks))
            loc = 1;
        end
        L = LaminsData(MD);
        L.params.steerable.sigma = 5;
        images = L.getImages;
        S = cell(numel(images),1);
        S2 = cell(numel(images),1);
        for ii=loc'
            try
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
            catch err
                if(isvalid(h))
                    close(h);
                end
                warning(['Failure ' MD.movieDataFileName_ ' for item ' num2str(ii)]);
                disp(err);
            end
       end
       file = [MD.outputDirectory_ filesep 'skeletons_20150519.mat'];
       save(file,'S','S2');
    end

end

