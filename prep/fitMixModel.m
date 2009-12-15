function [estimates] = fitMixModel(image, xyvec, sig, amps, b)
% fit mixture model to image with Jacobian to retrieve subpixel values for
% speckle positions, using info from cands as start point
% INPUT:    image
%           xyvec: (numspec x 2) vector containing startpoints for x,y,
%                   coordinates of speckles
%           sig:    sigma of Guassian, corresponding to (1/3)*PSF width
%           amps:   amplitudes (deltaI) vector
%           b:      background value
% OUTPUT: estimates
% UPDATED Dec 10/2009 Sylvain Berlemont

[xs,ys] = size(image);
numspec = size(xyvec, 1);

% IN ORDER TO FIT USING SPARSE MATRIX, FIRST IDENTIFY RELEVANT AREAS
% i.e., the ones SURROUNDING SPECKLES
% make mask of image with value 1 in a (3*sig radius) circle around
% speckles and value 0 elsewhere
% an additional calculation is made using 6*sig, to provide the templates
% for the superposition of single-speckle contributions to the the fit

%mask is copy of image with 1 for speckle vicinity; max projection of all
%single mat3's (calculated in loop below)
mask=zeros(xs, ys);

%masksin = sparse matrix with 1/0 values for coordinates; each speckle is
%one column
%masksin=sparse(xs*ys,numspec);

% radius for speckle spot calculation
rad = ceil(6*sig);
[miniImX, miniImY] = ndgrid(-rad:rad, -rad:rad);
miniDist = sqrt(miniImX.^2 + miniImY.^2);
%inside pixels
[XinPix, YinPix] = find(miniDist<=(3*sig));
XinPix = XinPix-(rad+1);
YinPix = YinPix-(rad+1);

% same for larger radius
[XinPixLar, YinPixLar] = find(miniDist<=(6*sig));
XinPixLar = XinPixLar-(rad+1);
YinPixLar = YinPixLar-(rad+1);

singleSpecLarPositions = cell(numspec,1);
GlobalCoordMat = cell(numspec,1);

hm = waitbar(0,'generating sparse mask');
for n=1:numspec
    SmalCoorX = xyvec(n, 1) + XinPix;
    SmalCoorY = xyvec(n, 2) + YinPix;
    % ind: positions inside the image
    ind = find((SmalCoorX>=1) & (SmalCoorY>=1) & (SmalCoorX<=xs) & (SmalCoorY<=ys));
    SmalCoorX = SmalCoorX(ind);
    SmalCoorY = SmalCoorY(ind);
    GlobalCoordMat{n} = sub2ind([xs ys], SmalCoorX, SmalCoorY);
    
    LarCoorX = xyvec(n,1) + XinPixLar;
    LarCoorY = xyvec(n,2) + YinPixLar;
    % ind: positions inside the image
    ind = find((LarCoorX>=1) & (LarCoorY>=1) & (LarCoorX<=xs) & (LarCoorY<=ys));
    LarCoorX = LarCoorX(ind);
    LarCoorY = LarCoorY(ind);
    singleSpecLarPositions{n} = [LarCoorX, LarCoorY, sub2ind([xs ys],LarCoorX, LarCoorY)]; %% ???
    
    waitbar(n / numspec, hm);
end

close(hm);

ind = cell2mat(GlobalCoordMat);
mask(ind) = 1;

% startpoint for fit
start_point = [xyvec amps];
lowbound=[xyvec-2 0.8*amps];
upbound=[xyvec+2 1.2*amps];
maxIter = 15;
options = optimset('TolFun',1e-4,'TolX',1e-4,'Display','off',...
    'MaxIter', maxIter,'Jacobian','on');

fh = waitbar(0,'fitting mixture model');
iteration = 1;

[estimates] = lsqnonlin(@Gauss2Dfun, start_point, lowbound, upbound, options);

% in the present implementation, the function Gauss2Dfun fits the x-
% and y-coordinates of speckles and the speckle amplitude.

% the returned function is the function difference F-image, NOT the squared
% error; the square minimization is done implicitly in the lsqnonlin
% function

    function [F,J] = Gauss2Dfun(params)
        x0 = params(:,1);
        y0 = params(:,2);
        A = params(:,3);
        ix = size(params, 1);
               
        %initialize fitimage
        Fitimage = zeros(size(image));

        %loop over number of speckles
        for z=1:ix
            %relevant 6*sig-distance positions of this speckle
            xyloc = singleSpecLarPositions{z}(:,1:2);
            iloc = singleSpecLarPositions{z}(:,3);
            zz=( (xyloc(:,1)-x0(z)).^2 + (xyloc(:,2)-y0(z)).^2 );
            % single speckle contribution
            FittedCurveSin = A(z) * exp( -(1/(2*(sig^2))) * zz);
            % contribution is added to Fitimage at appropriate position
            Fitimage(iloc) = Fitimage(iloc) + FittedCurveSin;
        end
        
        % multiply by mask (also pre-defined outside the actual fitting
        % function to save time) in order to take into account for the fit
        % only the relevant areas around speckles
        % and not try to fit unspecific background
        F = (Fitimage + b - image) .* mask;
        
        % now calculate the Jacobian matrix!!
        
        if (nargout>1)
            % fun returns vector of m=(xs*ys) components (corresponding to
            % image size), and parameter vector x has n=(numspec x 3) components,
            % namely the x and y coordinates for numspec speckles and their
            % amplitudes; thus, Jacobian is m x n matrix containing partial
            % derivatives for each differentiating parameter, fill up the
            % J-coordinates corresponding to speckle coordinates in the
            % image with the value one. The speckle coordinate positions
            % (3-sig circle) are stored in GlobalCoordMat calculated
            % outside the function.
            % The first parameter x=(x1,x2,x3,...xnumspec) corresponds to
            % the columns 1:numspec in the Jacobian, the second parameter:
            % y=(y1,y2,y3,....ynumspec) corresponds to columns
            % (numspec+1):end and same for third parameter, amplitude I0.

            % The following three variables are used to stored pair of
            % indices / value to build the Jacobian matrix J using:
            % J = sparse(e1, e2, e3, xs * ys, 3 * numspec);
            e1 = cell(numspec, 1);
            e2 = cell(numspec, 1);
            e3 = cell(numspec, 1);
            
            %now calculate actual values for Jacobian for each speckle
            for s=1:numspec
                %general procedure:
                %make left- and right-shifted image, subtract and divide by
                %two, according to
                %dx/dt=[(dx+Dt)-(dx-Dt)] / (2*Dt)
                
                %v non zero value position vector for speckle s
                vind = GlobalCoordMat{s};
                [xx,yy]=ind2sub([xs,ys],vind);
                
                %shift distance
                sh=0.3;
                
                %left shifted image of speckle s
                %not entire image is calculated, but just projection onto
                %nonzeros in v (due to structure of xx and yy above)
                zzleft=( (xx-x0(s)-sh).^2 + (yy-y0(s)).^2 );
                varleft = A(z) * exp( -(1/(2*(sig^2))) * ( zzleft) );
                timageleft=varleft-image(vind) + b;
                %right shifted image of speckle s
                zzright=( (xx-x0(s)+sh).^2 + (yy-y0(s)).^2 );
                varright = A(z) * exp( -(1/(2*(sig^2))) * ( zzright) );
                timageright=varright-image(vind) + b;
                %central difference image x (left-right)
                diff1=(-timageright+timageleft)/2;
                
                %up-shifted image of speckle s
                zztop=( (xx-x0(s)).^2 + (yy-y0(s)+sh).^2 );
                vartop = A(z) * exp( -(1/(2*(sig^2))) * ( zztop) );
                timagetop=vartop-image(vind) + b;
                %down-shifted image of speckle s
                zzbottom=( (xx-x0(s)).^2 + (yy-y0(s)-sh).^2 );
                varbottom = A(z) * exp( -(1/(2*(sig^2))) * ( zzbottom) );
                timagebottom=varbottom-image(vind) + b;
                %central difference image y (bottom-top)
                diff2=(-timagetop+timagebottom)/2;
                
                %differentiation df/da (amplitude)
                zzamp=( (xx-x0(s)).^2 + (yy-y0(s)).^2 );
                diff3=exp( -(1/(2*(sig^2))) * ( zzamp) );
                
                e1{s} = repmat(vind, 3, 1);
                e2{s} = vertcat(s * ones(size(vind)), ...
                    (s + numspec) * ones(size(vind)), ...
                    (s + 2 * numspec) * ones(size(vind)));
                e3{s} = vertcat(diff1, diff2, diff3);
            end
            
            J = sparse(vertcat(e1{:}), vertcat(e2{:}), vertcat(e3{:}), ...
                xs * ys, 3 * numspec);
            
        end % of if (nargout>1)
        iteration = iteration + 1;
        waitbar(iteration / maxIter, fh);
        
    end %of function F

close(fh);

end