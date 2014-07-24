function positions = loadPositions(projectNumber)
%========================================================
% POSITION DATA
%========================================================
% positions is a nTp-by-4-by-nTags array
% amplitudes are set to 1
% default movie size: 64 x 64 x 16 x 1 x nTP
% default pixel size: 0.05/0.05/0.2

switch projectNumber
    case {100}
        % Fixed spb, cen approaches vertically (horizontal=tag-txt overlp)
        % 13 frames
        % default size

        spbPosition = repmat([15,32,9,1],13,1);

        cenPosition = repmat([15,32,9,1],13,1);
        cenPosition(:,1) = cenPosition(:,1) + ...
            [20,20,20,12,8,7,6,5,4.5,4,3.5,3,2.5]';

        positions = cat(3,spbPosition,cenPosition);

    case {101}
        % Both SPB and CEN fixed.
        % 2 frames
        % default size
        spbPosition = repmat([15,32,9,1],2,1);

        cenPosition = repmat([32,32,9,1],2,1);

        positions = cat(3,spbPosition,cenPosition);
        
    case {102}
        % One fixed spot, six frames in a mini-movie (30x30x15x6)
        
        positions = repmat([15,15,8,1],[6,1]);
        
    case {103}
        
        % SPB fixed in center. CEN moves from -2hPSF to 2hPSF in
        % 0.1PSF-steps.
        % movieSize: standard
        
        % hPSF = 228nm
        fractPSF = 228/50;
        
        spbPosition = repmat([32,32,9,1],41,1);
        cenPosition = spbPosition;
        cenPosition(:,1) = cenPosition(:,1) +...
            fractPSF*(-2:0.1:2)';
        
        positions = cat(3,spbPosition,cenPosition);


    case {900}
        % one single spot exploring one eight of a voxel on a 0.1-voxel
        % grid in a mini-movie (20x20x10)
        [a,b,c]=ndgrid(10:0.1:10.5,10:0.1:10.5,5:0.1:5.5);

        positions = [a(:),b(:),c(:),ones(216,1)];

    case {901}
        % fixed SPB, CEN moves on iso-Rayleigh-surfaces
        % movieSize 35x35x24

        rayleigh = [3,2,1.5:-0.1:0.3];
        nR = 15;
        phi = [0,1/8,1/4] * pi;
        nP = 3;
        theta = [0:1/8:1/2] * pi;
        nT = 5;

        nm2pix = 1./[50,50,200];
        rayXY = 228.8;
        rayZ = 814.3;

        [spbPosition,cenPosition] = deal(repmat([10,10,7,1],225,1));
        i = 0;
        for t = 1:nT
            ct = cos(theta(t));
            st = sin(theta(t));
            rayT = 1/sqrt((ct/rayXY)^2 + (st/rayZ)^2);
            for p = 1:nP
                cp = cos(phi(p));
                sp = sin(phi(p));
                for r = 1:nR
                    % increment index
                    i = i + 1;
                    % Rayleigh-limit * factor
                    rayRT = rayT * rayleigh(r);
                    % angle-multiplicators
                    angles = [ct*cp, ct*sp, st];
                    cenPosition(i,1:3) = cenPosition(i,1:3) + ...
                        rayRT * angles .* nm2pix;
                end
            end
        end

        positions = cat(3,spbPosition,cenPosition);
        
    case 902
        
       % 16 spots, all far from each other. Positions randomly assume one
       % of the 216 points on the 1/8th voxel grid
       % 80x80x10x216 movie
       
       % remember random state
       randomState = rand('state');
       
       % set random state to 1
       rand('state',1);
       
       % find sub-voxel positions. Z will always be around 5
       [x,y,z]=ndgrid(0:0.1:0.5,0:0.1:0.5,5:0.1:5.5);
       
       % position matrix is 216 x 4 x 16
       positions = [zeros(216,3,16),ones(216,1,16)];
       
       % centers are at 10, 30, 50, and 70 pixels
       [u,v]=ndgrid(10:20:70);
       
       % loop through 16 spots and place
       for i=1:16
           positions(:,1,i) = x(randperm(216)') + u(i);
           positions(:,2,i) = y(randperm(216)') + v(i);
           positions(:,3,i) = z(randperm(216)');
       end
       
       % set random state back to where it was
       rand('state',randomState);
       

    otherwise
        disp(sprintf('project number %i is not defined',currentTest))
end
