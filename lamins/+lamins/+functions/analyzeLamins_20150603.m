function [ ] = analyzeLamins_20150603( MD )
%analyzeLamins_20150603 Batch analysis routine for June 3rd, 2015
% Adapted code to clear out large holes in LB1null mutants

    import lamins.functions.*;
    import lamins.classes.*;
    
    if(ischar(MD))
        try
            MD = MovieData.load(MD);
        catch err
            MD = MovieList.load(MD);
        end
    end

    if(isa(MD,'MovieList'))
        ML = MD;
        for mdi = 1:length(ML.movies_)
            MD = ML.getMovie(mdi);
            try
                analyzeLamins_20150603(MD);
            catch err
                warning(['Failure on : ' MD.movieDataFileName_]);
                disp(err);
            end
        end
        return;
    else
        R = CellReader(MD.getReader());
        C = R.toCell;
        % reorder so that we go in tzc order instead of ctz order
        C = shiftdim(C,2);
        means = cellfun(@(x) mean(x(:)),C(:,:));
        [pks,loc] = findpeaks(means(:));
        [tz,channel] = ind2sub(size(C),loc);
        tz = unique(tz);
        if(isempty(tz))
            tz = 1;
        end
        frames = zeros(length(tz)*R.getSizeC(),2);
        for c=1:R.getSizeC
            start = 1+(c-1)*length(tz);
            stop = start + length(tz) - 1;
            frames(start:stop,:) = [ones(size(tz))*c tz];
        end
        frames = sub2ind(R.getSize,frames(:,1),frames(:,2));
%         loc = 1;
        L = LaminsData(MD);
        L.params.steerable.sigma = 5;
        images = L.getImages;
        S = cell(numel(images),1);
        S2 = cell(numel(images),1);
        S3 = cell(numel(images),1);
        % note: probably had a bug in 20150519 analysis code
        for ii=frames'
            try
                I = images(ii);
                S{ii} = I.skeleton;
                S{ii}.cleanup
                score = lamins.functions.scoreEdges(S{ii},I.flattenIntensity);
                S{ii}.deleteEdges(score < 0);

                S2{ii} = S{ii}.copy;
                thresh = I.maskThresh(double(I));
                S2{ii}.auditEdges(double(I),[],thresh,thresh/2);
                
                % New on June 3rd, 2015
                % Audit using mask, proportion on flattened intensity, and
                % do another round of score optimization including zero
                S3{ii} = S2{ii}.copy;
                S3{ii}.auditEdgesByMask(I);
                S3{ii}.auditEdgesByThresholdedIntensity(I);
                score2 = lamins.functions.scoreEdges(S3{ii},I.flattenIntensity);
                S3{ii}.deleteEdges(score <= 0);

                h = figure;
                imshow(I);
                S{ii}.drawEdgesAsLines([],'r');
                S2{ii}.drawEdgesAsLines([],'g');
                S3{ii}.drawEdgesAsLines([],'b');
                title(MD.movieDataFileName_);
                text(50,50,[MD.movieDataFileName_ ' analyzed on 2015_06_03'],'Color','y','Interpreter','none')
                savefig(h,[MD.outputDirectory_ filesep 'skeleton_20150603_' num2str(ii) '.fig']);
                saveas(h, [MD.outputDirectory_ filesep 'skeleton_20150603_' num2str(ii) '.pdf']);
                saveas(h, [MD.outputDirectory_ filesep 'skeleton_20150603_' num2str(ii) '.png']);
                close(h);
            catch err
                if(isvalid(h))
                    close(h);
                end
                warning(['Failure ' MD.movieDataFileName_ ' for item ' num2str(ii)]);
                disp(err);
            end
       end
       file = [MD.outputDirectory_ filesep 'skeletons_20150603.mat'];
       save(file,'S','S2','S3');
    end

end

