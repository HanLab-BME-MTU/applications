function [ out ] = analyzeLaminsForProcessWithNLMS( MD , varargin)
%analyzeLaminsForProcess Batch analysis routine for lamin skeleton analysis
% Adapted code to clear out large holes in LB1null mutants
%
% MD is a MovieData object
%    OR
%    a MovieList object
%    a string for a MovieData or MovieList object

    % setup pointer back to this function
    fxn = @analyzeLaminsForProcessWithNLMS;
    out = [];

    import lamins.functions.*;
    import lamins.classes.*;
       
    % handle string input
    if(ischar(MD))
        try
            MD = MovieData.load(MD);
        catch
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
        ip = inputParser;
        ip.StructExpand = true;
        ip.KeepUnmatched = true;
        ip.addRequired('MD',@(x) isa(x,'MovieObject'));
        analyzeLaminsForProcessParameters(ip);
        ip.addParameter('process',[]);
        ip.parse(MD,varargin{:});
        
        in = ip.Results;
        
        if(isempty(in.channels))
            in.channels = 1:length(MD.channels_);
        else
            outOfRangeChannels = in.channels > length(MD.channels_);
            if(any(outOfRangeChannels))
                in.channels = in.channels(~outOfRangeChannels);
                warning('analyzeLaminsForProcess:outOfRangeChannels', ...
                    'channels parameter is out of range.');
            end
        end

        % setup date for saving files
        analysisDate = in.analysisDate;
    
        
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
        if(isempty(in.tz))
            tz = findFocalPlane(MD);
            tz = unique(tz);
            if(isempty(tz))
                tz = 1;
            end
        else
            tz = in.tz;
        end
        
%         [pks,loc] = findpeaks(means(:));
%         [tz,channel] = ind2sub(size(C),loc);

%         frames = zeros(length(tz)*length(in.channels),2);
%         for c=1:length(in.channels)
%             start = 1+(c-1)*length(tz);
%             stop = start + length(tz) - 1;
%             frames(start:stop,:) = [ones(size(tz))*in.channels(c) tz];
%         end
%         frames = sub2ind(R.getSize,frames(:,1),frames(:,2));
        
%         loc = 1;
        L = LaminsData(MD);
        L.params.steerable.sigma = in.steerable_sigma;
        images = L.getImages;
%         S = cell(numel(images),1);
%         S2 = cell(numel(images),1);
%         S3 = cell(numel(images),1);
        % note: probably had a bug in 20150519 analysis code
        for c = in.channels
            S = cell(max(tz),1);
            S2 = cell(max(tz),1);
            S3 = cell(max(tz),1);
            S4 = cell(max(tz),1);
            if(in.clearOutputDir)
                mkClrDir(in.outFilePaths{c});
            else
                if(~exist(in.outFilePaths{c},'dir'))
                    mkdir(in.outFilePaths{c});
                end
            end
            for tz_i = tz(:)'
                try
                    I = images(c,tz_i);
                    [S4{tz_i},S3{tz_i},S2{tz_i},S{tz_i}] = analyzeLaminsImage(I);
                    if(in.plot)
                        plotFigures(I, ...
                            S{tz_i},S2{tz_i},S3{tz_i},S4{tz_i},...
                            MD.movieDataFileName_, ...
                            MD.channels_(c).name_, ...
                            in.outFilePaths{c}, ...
                            analysisDate, ...
                            tz_i);
                    end
                catch err
                    % Close figure if there is a failure
                    if(exist('h','var') && isvalid(h))
                        close(h);
                    end
                    warning(['Failure ' MD.movieDataFileName_ ' for item ' num2str(tz_i) ' Channel ' num2str(c)]);
                    disp(getReport(err));
                end
            end
            
            file = [in.outFilePaths{c} filesep 'skeletons_' analysisDate '.mat'];
            if(exist(file,'file'))
                % Merge new results into old data
                data = load(file);
                data.S(tz) = S(tz);
                data.S2(tz) = S2(tz);
                data.S3(tz) = S3(tz);
                data.S4(tz) = S4(tz);
                data.tz = unique([data.tz(:) ; tz(:) ]);
                save(file,'-struct','data')
            else
                save(file,'S','S2','S3','S4','tz');
            end
        end
    end
%         for ii=frames'
%             try
%                 I = images(ii);
%                 [S3{ii},S2{ii},S{ii}] = analyzeLaminsImage(I);
%                 plotFigures(I, ...
%                     S{ii},S2{ii},S3{ii},...
%                     MD.movieDataFileName_, ...
%                     MD.channels_(I.coordinates{1}).name_, ...
%                     MD.outputDirectory_, ...
%                     analysisDate, ...
%                     ii);
%                 % Initial skeleton with basic cleanup
%                 S{ii} = I.skeleton;
%                 S{ii}.cleanup
%                 score = lamins.functions.scoreEdges(S{ii},I.flattenIntensity);
%                 S{ii}.deleteEdges(score < 0);
% 
%                 % Use intensity variation along the edge
%                 S2{ii} = S{ii}.copy;
%                 thresh = I.maskThresh(double(I));
%                 S2{ii}.auditEdges(double(I),[],thresh,thresh/2);
%                 
%                 % New on June 3rd, 2015
%                 % Audit using mask, proportion on flattened intensity, and
%                 % do another round of score optimization including zero
%                 S3{ii} = S2{ii}.copy;
%                 S3{ii}.auditEdgesByMask(I);
%                 S3{ii}.auditEdgesByThresholdedIntensity(I);
%                 score2 = lamins.functions.scoreEdges(S3{ii},I.flattenIntensity);
%                 S3{ii}.deleteEdges(score2 <= 0);
% 
%                 % Figure 1: Illustrate progression of auditing procedures
%                 h = figure;
%                 imshow(im2uint8(I.adjusted));
%                 S{ii}.drawEdgesAsLines([],'r');
%                 S2{ii}.drawEdgesAsLines([],'g');
%                 S3{ii}.drawEdgesAsLines([],'b');
%                 title(MD.movieDataFileName_);
%                 text(50,50,[MD.movieDataFileName_ ' analyzed on ' analysisDate],'Color','y','Interpreter','none')
%                 savefig(h,[MD.outputDirectory_ filesep 'skeleton_' analysisDate '_' num2str(ii) '.fig']);
%                 saveas(h, [MD.outputDirectory_ filesep 'skeleton_' analysisDate '_' num2str(ii) '.pdf']);
%                 saveas(h, [MD.outputDirectory_ filesep 'skeleton_' analysisDate '_' num2str(ii) '.png']);
%                 close(h);
%                 
%                 % Figure 2: Show final skeleton on reverse complement
%                 h = figure;
%                 imshow(im2uint8(I.adjustedComplement));
%                 S3{ii}.drawEdgesAsLines([],'m');
%                 title(MD.movieDataFileName_);
%                 text(50,50,[MD.movieDataFileName_ ' analyzed on ' analysisDate],'Color','b','Interpreter','none')
%                 savefig(h,[MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(ii) '.fig']);
%                 saveas(h, [MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(ii) '.pdf']);
%                 saveas(h, [MD.outputDirectory_ filesep 'one_skeleton_' analysisDate '_' num2str(ii) '.png']);
%                 close(h);
%                 
%                 % Figure 3: Show image by itself
%                 h = figure;
%                 imshow(im2uint8(I.adjusted));
%                 title(MD.movieDataFileName_);
%                 text(50,50,[MD.movieDataFileName_ ' viewed on ' analysisDate],'Color','y','Interpreter','none');
%                 text(50,80,['Channel ' MD.channels_(I.coordinates{1}).name_],'Color','y','Interpreter','none');
%                 text(50,110,['Z-slice ' num2str(I.coordinates{3}) ', Frame ' num2str(I.coordinates{2})],'Color','y','Interpreter','none');
%                 saveas(h, [MD.outputDirectory_ filesep 'image_' analysisDate '_' num2str(ii) '.png']);
%                 close(h);
%                 
%                 % collect properties of faces for Figures 4 and beyond
%                 rp = regionprops(S3{ii}.faces,'Area','Eccentricity','Perimeter');
%                 
%                 % Figure 4: Histogram of face area
%                 h = figure;
%                 histogram([rp.Area]);
%                 ylabel('Number of Faces');
%                 xlabel('Area (\mum^2)');
%                 savefig(h,[MD.outputDirectory_ filesep 'Face_Area_Histogram_' analysisDate '_' num2str(ii) '.fig']);
%                 saveas(h,[MD.outputDirectory_ filesep 'Face_Area_Histogram_' analysisDate '_' num2str(ii) '.png']);
%                 close(h);
%                 
%                 % Figure 5: Histogram of face eccentricity
%                 h = figure;
%                 histogram([rp.Eccentricity]);
%                 ylabel('Number of Faces');
%                 xlabel('Eccentricity');
%                 savefig(h,[MD.outputDirectory_ filesep 'Face_Eccentricity_Histogram_' analysisDate '_' num2str(ii) '.fig']);
%                 saveas(h,[MD.outputDirectory_ filesep 'Face_Eccentricity_Histogram_' analysisDate '_' num2str(ii) '.png']);
%                 close(h);
%                 
%                 % Figure 6: Face circularity
%                 h = figure;
%                 histogram(4*pi*[rp.Area]./[rp.Perimeter].^2);
%                 ylabel('Number of Faces');
%                 xlabel('Face Circularity, 1 = circle, 0 = line');
%                 savefig(h,[MD.outputDirectory_ filesep 'Face_Circularity_Histogram_' analysisDate '_' num2str(ii) '.fig']);
%                 saveas(h,[MD.outputDirectory_ filesep 'Face_Circularity_Histogram_' analysisDate '_' num2str(ii) '.png']);
%                 close(h);
%             catch err
%                 % Close figure if there is a failure
%                 if(exist('h','var') && isvalid(h))
%                     close(h);
%                 end
%                 warning(['Failure ' MD.movieDataFileName_ ' for item ' num2str(ii)]);
%                 disp(getReport(err));
%             end
%        end
%        file = [MD.outputDirectory_ filesep 'skeletons_' analysisDate '.mat'];
%        save(file,'S','S2','S3','tz');
%     end
    
    function out = errorHandler(errStruct, MD)
        out = [];
        warning(['Error in movie number ' num2str(errStruct.index) ' MD.movieDataFileName_']);
        warning(errStruct.message);
    end

end
function [S4,S3,S2,S] = analyzeLaminsImage(laminsImage)
%     I = images(ii);
    I = laminsImage;
    % Initial skeleton with basic cleanup
%     S = I.skeleton;
    F8 = OrientationSpaceFilter(0.08,0.08*0.8,8);
    R8 = F8*double(laminsImage.image);
    R3 = R8.getResponseAtOrderFT(3);
    maxima = R8.getRidgeOrientationLocalMaxima;
    nlms = nonLocalMaximaSuppressionPrecise(real(R3.a),maxima);
    nlms_mip = nanmax(nlms,[],3);
    S = lamins.classes.Skeleton(nlms_mip);
    S.cleanup
    score = lamins.functions.scoreEdges(S,I.flattenIntensity);
    S.deleteEdges(score < 0);

    % Use intensity variation along the edge
    S2 = S.copy;
    thresh = I.maskThresh(double(I));
    S2.auditEdges(double(I),[],thresh,thresh/2);

    % New on June 3rd, 2015
    % Audit using mask, proportion on flattened intensity, and
    % do another round of score optimization including zero
    % This was the version published with the Dec 2015 MBoC Paper
    S3 = S2.copy;
    S3.auditEdgesByMask(I);
    S3.auditEdgesByThresholdedIntensity(I);
    score2 = lamins.functions.scoreEdges(S3,I.flattenIntensity);
    S3.deleteEdges(score2 <= 0);
    
    % New on February 1st, 2016
    S4 = S3.copy;
    % Some vertices were not deleted resulting in shorter edges
    S4.mergeEdgesAtObsoleteVertices;
end
function plotFigures(laminsImage,S,S2,S3,S4,fileName,channelName,outputDirectory,analysisDate,num)

    if(isnumeric(num))
        num = num2str(num);
    end
    
    I = laminsImage;
    
    if(isa(fileName,'MovieData'))
        MD = fileName;
        fileName = MD.movieDataFileName_;
        outputDirectory = MD.outputDirectory_;
        channelName = MD.channels_(I.coordinates{1}).name_;
    end

    % Figure 1: Illustrate progression of auditing procedures
    h = figure;
    imshow(im2uint8(I.adjusted));
    S.drawEdgesAsLines([],'r');
    S2.drawEdgesAsLines([],'g');
    S4.drawEdgesAsLines([],'b');
    title(fileName);
    text(50,50,[fileName ' analyzed on ' analysisDate],'Color','y','Interpreter','none')
    savefig(h,[outputDirectory filesep 'skeleton_' analysisDate '_' num '.fig']);
    saveas(h, [outputDirectory filesep 'skeleton_' analysisDate '_' num '.pdf']);
    saveas(h, [outputDirectory filesep 'skeleton_' analysisDate '_' num '.png']);
    close(h);

    % Figure 2: Show final skeleton on reverse complement
    h = figure;
    imshow(im2uint8(I.adjustedComplement));
    S4.drawEdgesAsLines([],'m');
    title(fileName);
    text(50,50,[fileName ' analyzed on ' analysisDate],'Color','b','Interpreter','none')
    savefig(h,[outputDirectory filesep 'one_skeleton_' analysisDate '_' num '.fig']);
    saveas(h, [outputDirectory filesep 'one_skeleton_' analysisDate '_' num '.pdf']);
    saveas(h, [outputDirectory filesep 'one_skeleton_' analysisDate '_' num '.png']);
    close(h);

    % Figure 3: Show image by itself
    h = figure;
    imshow(im2uint8(I.adjusted));
    title(fileName);
    text(50,50,[fileName ' viewed on ' analysisDate],'Color','y','Interpreter','none');
    text(50,80,['Channel ' channelName],'Color','y','Interpreter','none');
    text(50,110,['Z-slice ' num2str(I.coordinates{3}) ', Frame ' num2str(I.coordinates{2})],'Color','y','Interpreter','none');
    saveas(h, [outputDirectory filesep 'image_' analysisDate '_' num '.png']);
    close(h);

    % collect properties of faces for Figures 4 and beyond
    rp = regionprops(S4.faces,'Area','Eccentricity','Perimeter');

    % Figure 4: Histogram of face area
    h = figure;
    histogram([rp.Area]);
    ylabel('Number of Faces');
    xlabel('Area (\mum^2)');
    savefig(h,[outputDirectory filesep 'Face_Area_Histogram_' analysisDate '_' num '.fig']);
    saveas(h,[outputDirectory filesep 'Face_Area_Histogram_' analysisDate '_' num '.png']);
    close(h);

    % Figure 5: Histogram of face eccentricity
    h = figure;
    histogram([rp.Eccentricity]);
    ylabel('Number of Faces');
    xlabel('Eccentricity');
    savefig(h,[outputDirectory filesep 'Face_Eccentricity_Histogram_' analysisDate '_' num '.fig']);
    saveas(h,[outputDirectory filesep 'Face_Eccentricity_Histogram_' analysisDate '_' num '.png']);
    close(h);

    % Figure 6: Face circularity
    h = figure;
    histogram(4*pi*[rp.Area]./[rp.Perimeter].^2);
    ylabel('Number of Faces');
    xlabel('Face Circularity, 1 = circle, 0 = line');
    savefig(h,[outputDirectory filesep 'Face_Circularity_Histogram_' analysisDate '_' num '.fig']);
    saveas(h,[outputDirectory filesep 'Face_Circularity_Histogram_' analysisDate '_' num '.png']);
    close(h);
end

