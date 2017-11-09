function [ out ] = analyzeLamins_20150610( MD )
%analyzeLamins_20150603 Batch analysis routine for June 3rd, 2015
% Adapted code to clear out large holes in LB1null mutants
%
% MD is a MovieData object
%    OR
%    a MovieList object
%    a string for a MovieData or MovieList object

    out = [];

    import lamins.functions.*;
    import lamins.classes.*;
    
    % setup date for saving files
    analysisDate = '2015_06_10';
    % setup pointer back to this function
    fxn = @analyzeLamins_20150610;
    
    % handle string input
    if(ischar(MD))
        try
            MD = MovieData.load(MD);
        catch err
            MD = MovieList.load(MD);
        end
    end

    if(isa(MD,'MovieList'))
        % if MovieList, then iterate fxn over all the movies
        ML = MD;
        movies = distributed(ML.movies_);
        cellfun(fxn,movies,'UniformOutput',false,'ErrorHandler',@errorHandler);
        
%         for mdi = 1:length(ML.movies_)
%             MD = ML.getMovie(mdi);
%             try
%                 fxn(MD);
%             catch err
%                 warning(['Failure on : ' MD.movieDataFileName_]);
%                 disp(getReport(err));
%             end
%         end
        return;
    else
        % Assume this is a MovieData object and analyze
        disp(MD.movieDataPath_);
        R = CellReader(MD.getReader());
%         C = R.toCell;
        % reorder so that we go in tzc order instead of ctz order
%         C = shiftdim(C,2);
        
%         % Perform autofocus based on edge detection values
%         h = -fspecial('log',[5 5],0.5);
% %         means = cellfun(@(x) mean(x(:)),C(:,:));
%         focusMetric = cellfun(@(x) sum(joinColumns(imfilter(double(x),h))),C(:,:));
%         [~,tz] = max(focusMetric);
%         tz = tz';
        tz = findFocalPlane(MD);
        
%         [pks,loc] = findpeaks(means(:));
%         [tz,channel] = ind2sub(size(C),loc);
%         tz = unique(tz);
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
                % Initial skeleton with basic cleanup
                S{ii} = I.skeleton;
                S{ii}.cleanup
                score = lamins.functions.scoreEdges(S{ii},I.flattenIntensity);
                S{ii}.deleteEdges(score < 0);

                % Use intensity variation along the edge
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
                S3{ii}.deleteEdges(score2 <= 0);

                % Figure 1: Illustrate progression of auditing procedures
                h = figure;
                imshow(im2uint8(I.adjusted));
                S{ii}.drawEdgesAsLines([],'r');
                S2{ii}.drawEdgesAsLines([],'g');
                S3{ii}.drawEdgesAsLines([],'b');
                title(MD.movieDataFileName_);
                text(50,50,[MD.movieDataFileName_ ' analyzed on ' analysisDate],'Color','y','Interpreter','none')
                savefig(h,[MD.outputDirectory_ filesep 'skeleton_' analysisDate '_' num2str(ii) '.fig']);
                saveas(h, [MD.outputDirectory_ filesep 'skeleton_' analysisDate '_' num2str(ii) '.pdf']);
                saveas(h, [MD.outputDirectory_ filesep 'skeleton_' analysisDate '_' num2str(ii) '.png']);
                close(h);
                
                % Figure 2: Show final skeleton on reverse complement
                h = figure;
                imshow(im2uint8(I.adjustedComplement));
                S3{ii}.drawEdgesAsLines([],'m');
                title(MD.movieDataFileName_);
                text(50,50,[MD.movieDataFileName_ ' analyzed on ' analysisDate],'Color','b','Interpreter','none')
                savefig(h,[MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(ii) '.fig']);
                saveas(h, [MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(ii) '.pdf']);
                saveas(h, [MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(ii) '.png']);
                close(h);
                
                % Figure 3: Show image by itself
                h = figure;
                imshow(im2uint8(I.adjusted));
                title(MD.movieDataFileName_);
                text(50,50,[MD.movieDataFileName_ ' viewed on ' analysisDate],'Color','y','Interpreter','none');
                text(50,80,['Channel ' MD.channels_(I.coordinates{1}).name_],'Color','y','Interpreter','none');
                text(50,110,['Z-slice ' num2str(I.coordinates{3}) ', Frame ' num2str(I.coordinates{2})],'Color','y','Interpreter','none');
                saveas(h, [MD.outputDirectory_ filesep 'image_' analysisDate '_' num2str(ii) '.png']);
                close(h);
                
                % collect properties of faces for Figures 4 and beyond
                rp = regionprops(S3{ii}.faces,'Area','Eccentricity','Perimeter');
                
                % Figure 4: Histogram of face area
                h = figure;
                histogram([rp.Area]);
                ylabel('Number of Faces');
                xlabel('Area (\mum^2)');
                savefig(h,[MD.outputDirectory_ filesep 'Face_Area_Histogram_' analysisDate '_' num2str(ii) '.fig']);
                saveas(h,[MD.outputDirectory_ filesep 'Face_Area_Histogram_' analysisDate '_' num2str(ii) '.png']);
                close(h);
                
                % Figure 5: Histogram of face eccentricity
                h = figure;
                histogram([rp.Eccentricity]);
                ylabel('Number of Faces');
                xlabel('Eccentricity');
                savefig(h,[MD.outputDirectory_ filesep 'Face_Eccentricity_Histogram_' analysisDate '_' num2str(ii) '.fig']);
                saveas(h,[MD.outputDirectory_ filesep 'Face_Eccentricity_Histogram_' analysisDate '_' num2str(ii) '.png']);
                close(h);
                
                % Figure 6: Face circularity
                h = figure;
                histogram(4*pi*[rp.Area]./[rp.Perimeter].^2);
                ylabel('Number of Faces');
                xlabel('Face Circularity, 1 = circle, 0 = line');
                savefig(h,[MD.outputDirectory_ filesep 'Face_Circularity_Histogram_' analysisDate '_' num2str(ii) '.fig']);
                saveas(h,[MD.outputDirectory_ filesep 'Face_Circularity_Histogram_' analysisDate '_' num2str(ii) '.png']);
                close(h);
            catch err
                % Close figure if there is a failure
                if(exist('h','var') && isvalid(h))
                    close(h);
                end
                warning(['Failure ' MD.movieDataFileName_ ' for item ' num2str(ii)]);
                disp(getReport(err));
            end
       end
       file = [MD.outputDirectory_ filesep 'skeletons_' analysisDate '.mat'];
       save(file,'S','S2','S3','tz');
    end
    
    function out = errorHandler(errStruct, MD)
        out = [];
        warning(['Error in movie number ' num2str(errStruct.index) ' MD.movieDataFileName_']);
        warning(errStruct.message);
    end

end

