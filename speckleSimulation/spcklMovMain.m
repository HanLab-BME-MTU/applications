function [frmMeanI,fluorData,params]=spcklMovMain(params)
% here we generate a movie from user-defined model of speckle motion

% NOTE: this function gets executed by spcklMovParams.  (The
% function names are somewhat deceptive - the params function is the one to
% run, not this one.) 


% subfunction to set set up directories and generate initial fluorophores
[params,not2keep,fluorData,warningstate]=initializeEverything(params);

% ===== MODEL SPECIFIC SECTION =====
% in this section you can change the STATE distributions

switch params.nModel
    case 1 % Stationary network, random poly/depoly

        % if state is 0, the fluorophore is unbound (in G-actin state)and does
        % not contribute to signal; if state is 1, the fluorophore is bound to
        % the actin network (in F-actin state) and does contribute.
        randoms=rand(size(fluorData.py(:,1)));

        % here we initialize half the fluorophores to be unbound, half bound
        fluorData.state(randoms<0.5,1) = 0;
        fluorData.state(randoms>=0.5,1)= 1;

        % calcuate all future positions, state changes, etc.
        [params,fluorData,not2keep]=spcklMovModel1(params,fluorData,not2keep);

    case 2 % Converging network, either without depoly or depoly so I is constant

        % in this model we start with every labeled fluorophore bound
        fluorData.state(:)=1;

        % calcuate all future positions, state changes, etc.
        [params,fluorData,not2keep]=spcklMovModel2(params,fluorData,not2keep);

    case 3 % single-direction flow actin (one color) or actin/adhesion (two colors)
        
        randoms=rand(size(fluorData.py(:,1)));
        
        if params.protein == 1 % model actin
            % here we initialize half the fluorophores to be unbound, half
            % bound
            fluorData.state(randoms<0.5,1) = 0;
            fluorData.state(randoms>=0.5,1)= 2;

        elseif params.protein ==2 % model adhesion
            % here we initialize half the fluorophore states depending on
            % user-selected state distributions
            switch params.stateDist
                case 1 % all unbound
                    fluorData.state(:,1)=0;
                case 2 % all bound to substrate
                    fluorData.state(:,1)=1;
                case 3 % all bound to actin
                    fluorData.state(:,1)=2;
                case 4 % all bound to both substrate and actin
                    fluorData.state(:,1)=3;
                case 5 % random even distribution of monomers in the four states
                    fluorData.state(randoms<0.25,1)=0;
                    fluorData.state(randoms>=0.25 & randoms<0.50,1)=1;
                    fluorData.state(randoms>=0.50 & randoms<0.75,1)=2;
                    fluorData.state(randoms>=0.75,1)=3;
                case 6 % don't allow binding to both substrate and actin simultaneously
                    fluorData.state(randoms<0.33,1)=0;
                    fluorData.state(randoms>=0.33 & randoms<0.66,1)=1;
                    fluorData.state(randoms>=0.66,1)=2;
                otherwise
            end
        end
        
        % calcuate all future positions, state changes, etc.
        [params,fluorData,not2keep]=spcklMovModel3(params,fluorData,not2keep);

    case 4 % Antiparallel microtubules in bundles

        fluorData.state(:)=1;
        
        % calcuate all future positions
        [params,fluorData,not2keep]=spcklMovModel4(params,fluorData,not2keep);

    otherwise
end
% ===== END OF MODEL SPECIFIC SECTION =====


% FIGURE OUT MORE STUFF FOR fluorData

% .pixY and .pixX are the yx-indices of the pixel in which each fluorophore
% falls at each time point, including border pixels.FluorData.pixY=ceil(fluorData.py./params.pixNM);
fluorData.pixY=ceil(fluorData.py./params.pixNM);
fluorData.pixX=ceil(fluorData.px./params.pixNM);
% to make sure each falls inside a pixel, round
fluorData.pixY(fluorData.pixY>not2keep.imgLwB)=not2keep.imgLwB;
fluorData.pixX(fluorData.pixX>not2keep.imgWwB)=not2keep.imgWwB;

% .inView records whether the fluorophore is inside the field of view (1)
% or on the border pixels (0). this will be useful later for determining
% accurate turnover info, for example.
fluorData.inView(fluorData.pixY>=params.border.top+1 & ...
    fluorData.pixY<=params.border.top+params.imL & ...
    fluorData.pixX>=params.border.left+1 & ...
    fluorData.pixX<=params.border.left+params.imW)=1;


% WRITE THE IMAGE
[frmMeanI]=spcklMovWrtIm(params,fluorData,not2keep);

% SAVE THE IMPORTANT STUFF IN THE 'MovieInfo' DIRECTORY
save([params.outputDirMovieInfo filesep 'frmMeanI'],'frmMeanI');
save([params.outputDirMovieInfo filesep 'parameters'],'params');
save([params.outputDirMovieInfo filesep 'fluorStateAndLoc'],'fluorData');
save([params.outputDirMovieInfo filesep 'not2keep'],'not2keep'); % oh, why not?

% RETURN THE WARNING SETTINGS TO THEIR ORIGINAL STATE
warning(warningstate);
end



% SUBFUNCTIONS ================================================
function [params,not2keep,fluorData,warningstate]=initializeEverything(params)
% this subfunction includes everything that should be done before the main
% iterations are done according to model type.  it exists here so that the part
% of the function needing to be changed by the user will be visible at the top

% SAVE THE ORIGINAL WARNING SETTINGS AND TURN ANNOYING ONES OFF
warningstate=warning;
warning off MATLAB:intConvertNonIntVal

% SET UP THE DIRECTORIES ACCORDING TO THE MODEL TYPE
% this is a subfunction
[params]=spcklMovDirectories(params);


% GET TOP, BOTTOM, LEFT, AND RIGHT BORDERS FOR FRAMES
[params,not2keep]=spcklMovCalcBorder(params);


% INITIALIZE FLUOROPHORES AT T=0
[params,py1,px1] = spckleMovGenFluor(params,not2keep);


% CALCULATE STORAGE INFO
% num time points in a frame
not2keep.nTmPtsPerFm=params.nSecPerFrame/params.dT;

% num time points stored per frame
not2keep.nTmPtsStoredPerFm=round(params.intTime/params.dT);

% num time points in whole movie
not2keep.nTmPtsInMovie=not2keep.nTmPtsPerFm*(params.nFrames-1)+not2keep.nTmPtsStoredPerFm;

% num time points stored over whole movie
not2keep.nTmPtsStoredInMovie=not2keep.nTmPtsStoredPerFm*params.nFrames;

% points to store - vector of frame numbers where data should be stored
not2keep.pts2Store=repmat(1:not2keep.nTmPtsStoredPerFm,[params.nFrames 1])+...
    repmat(not2keep.nTmPtsPerFm*[0:params.nFrames-1]',[1,not2keep.nTmPtsStoredPerFm]);
not2keep.pts2Store=not2keep.pts2Store';
not2keep.pts2Store=not2keep.pts2Store(:);


% INITIALIZE "fluorData" STRUCTURE to store info about every fluorophore at
% all time points within integration time each frame
%   .state:     meaning varies by model, but used to determine speed, if visible or
%               diffuse, etc.
%   .py(x):     y(x)-coord (nm) for each fluorophore at each stored time point
%   .pixY(X):   y(x)-coord (pix) for each fluorophore at each stored time
%               point
%   .inView:    keeps track of whether each fluorophore is in the
%               field of view (1) or on a boundary pixel (0)
temp=zeros(length(px1),not2keep.nTmPtsStoredInMovie);
fluorData.state=temp;
fluorData.py=temp;
fluorData.px=temp;
fluorData.pixY=temp;
fluorData.pixX=temp;
fluorData.inView=temp;

% put coordinates from genFluor into the first time point
fluorData.py(:,1)=py1;
fluorData.px(:,1)=px1;

end % end of initializeEverything subfunction

% -------------------------------------------------------------------------

function [params]=spcklMovDirectories(params)
% this function sets up the directories for each model appropriately

switch params.nModel

    case {1, 2, 4}

        % /top-level project folder (user-selected)
        %   /analysis
        %       /edge/cell_mask
        %   /images/tifs
        %   /movieInfo
        %
        params.projDir=uigetdir(pwd,'Create or choose project directory'); % top level
        if isdir(params.projDir)
            rmdir(params.projDir,'s')
        end
        mkdir(params.projDir);

        params.anDir=[params.projDir filesep 'analysis' filesep 'edge'  filesep 'cell_mask'];
        mkdir(params.anDir);

        params.outputDirTifs=[params.projDir filesep 'images' filesep 'tifs'];
        mkdir(params.outputDirTifs);

        params.outputDirMovieInfo=[params.projDir filesep 'movieInfo'];
        mkdir(params.outputDirMovieInfo);

    case 3
        
        switch params.protein
            case 1 % actin

                % /top-level project folder (user-selected)
                %   /red_actin
                %       /analysis
                %           /edge
                %               /cell_mask
                %       /images
                %           /tifs
                %       /movieInfo

                params.projDir=uigetdir(pwd,'Create or choose top-level project directory'); % top level

                params.actinDir=[params.projDir filesep 'red_actin'];
                if isdir(params.actinDir)
                    rmdir(params.actinDir,'s')
                end
                mkdir(params.actinDir);

                params.anDir=[params.actinDir filesep 'analysis' filesep 'edge'  filesep 'cell_mask'];
                mkdir(params.anDir);

                params.outputDirTifs=[params.actinDir filesep 'images' filesep 'tifs'];
                mkdir(params.outputDirTifs);

                params.outputDirMovieInfo=[params.actinDir filesep 'movieInfo'];
                mkdir(params.outputDirMovieInfo);

            case 2 % adhesion

                % /top-level project folder (user-selected)
                %   /green_adhesion
                %       /analysis
                %           /edge/cell_mask
                %       /images/tifs
                %       /movieInfo
                %   /red_actin
                %       /analysis
                %           /edge/cell_mask
                %       /images/tifs
                %       /movieInfo

                params.projDir=uigetdir(pwd,'Choose top-level project directory (above /red-actin)'); % top level
                params.actinDir=[params.projDir filesep 'red_actin' filesep 'images' filesep 'tifs'];
                if ~isdir(params.actinDir)
                    error('Please run model 3 for actin first')
                end

                params.adhesDir=[params.projDir filesep 'green_adhesion'];
                if isdir(params.adhesDir)
                    rmdir(params.adhesDir,'s')
                end
                mkdir(params.adhesDir);

                params.anDir=[params.adhesDir filesep 'analysis' filesep 'edge'  filesep 'cell_mask'];
                mkdir(params.anDir);

                params.outputDirTifs=[params.adhesDir filesep 'images' filesep 'tifs'];
                mkdir(params.outputDirTifs);

                params.outputDirMovieInfo=[params.adhesDir filesep 'movie_info'];
                mkdir(params.outputDirMovieInfo);

            otherwise
                error('User input for model 3 should be 1 for actin or 2 for adhesion. Must run actin before adhesion if doing dual-channel.')
        end
    otherwise
end % end switch between model numbers for directory set-up

end % end of spcklMovDirectories subfunction





