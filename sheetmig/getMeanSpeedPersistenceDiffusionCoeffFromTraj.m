function [meanSpeedCell, persistenceCell, diffCoeffCell] = getMeanSpeedPersistenceDiffusionCoeffFromTraj(curFolder)
%% get the file list
fL = dir(curFolder);
fL = fL(3:end); % excluding . and ..
pixSize = 0.645; % um/pixel
numFiles = numel(fL);
%% preparing output arguments
meanSpeedCell = cell(numFiles,1);
persistenceCell = cell(numFiles,1);
diffCoeffCell = cell(numFiles,1);
%% going through each file
for ii=1:numFiles
    % Load your data
    dataOrg = readtable([fL(ii).folder filesep fL(ii).name]); 

    % Assuming data is in the format: [x1, y1, x2, y2, ...] for each time point
    
    % Number of cells
    num_cells = max(dataOrg.TrackN_);
    
    % Time step (you should replace this with your actual time step)
    dt = 10; % min 
    
    % Initialize arrays to store results
    mean_speeds = zeros(num_cells, 1);
    persistences = zeros(num_cells, 1);
    diffusion_coeffs = zeros(num_cells, 1);
    
    for cell_idx = 1:num_cells
        indCurCell = dataOrg.TrackN_ == cell_idx;
        % Extract x and y positions for the current cell
        x = dataOrg.X(indCurCell)*pixSize;
        y = dataOrg.Y(indCurCell)*pixSize;
        
        % Calculate displacements
        dx = diff(x);
        dy = diff(y);
        displacements = sqrt(dx.^2 + dy.^2);
        
        % Mean speed
        mean_speed = mean(displacements) / dt;
        mean_speeds(cell_idx) = mean_speed;
        
        % Persistence
        dot_products = dx(1:end-1).*dx(2:end) + dy(1:end-1).*dy(2:end);
        norms = displacements(1:end-1) .* displacements(2:end);
        cos_theta = dot_products ./ norms;
        persistence = mean(cos_theta(~isnan(cos_theta)));
        persistences(cell_idx) = persistence;
        
        % Diffusion coefficient (MSD/4t)
        MSD = mean(displacements.^2);
        diffusion_coeff = MSD / (4 * dt);
        diffusion_coeffs(cell_idx) = diffusion_coeff;
    end
    
    % Display results
    disp(['Mean Speeds (file ' num2str(ii) '):']);
    disp(mean_speeds);
    
    disp('Persistences:');
    disp(persistences);
    
    disp('Diffusion Coefficients:');
    disp(diffusion_coeffs);
    
    % Save results
    meanSpeedCell{ii} = mean_speeds;
    persistenceCell{ii} = persistences;
    diffCoeffCell{ii} = diffusion_coeffs;
end