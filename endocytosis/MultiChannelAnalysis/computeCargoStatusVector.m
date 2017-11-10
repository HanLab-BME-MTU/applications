function data = computeCargoStatusVector(data, channel, reference, parname_c0, parname_c1)
%
%
% Francois Aguet, Jan 2010

if (nargin <= 2)
    reference = 0;
end;
if (nargin>3 && ischar(parname_c0))
    paramName_c0 = [parname_c0 '.mat'];
    paramNameRef_c0 = [parname_c0 '_Ref.mat'];
else
    paramName_c0 = 'parameterMat.mat';
    paramNameRef_c0 = 'parameterMat_Ref.mat';
end;
if (nargin>4 && ischar(parname_c1))
    paramName_c1 = [parname_c1 '.mat'];
    paramNameRef_c1 = [parname_c1 '_Ref.mat'];
else
    paramName_c1 = 'SA_red_ParameterMat_.mat';
    paramNameRef_c1 = 'SA_red_ParameterMat__Ref.mat';
end;

nMovies = length(data);
fileName_c0 = cell(1:nMovies);
fileNameRef_c0 = cell(1:nMovies);
fileName_c1 = cell(1:nMovies);
fileNameRef_c1 = cell(1:nMovies);

%================================================================
% 1. Set parameter file names
%================================================================
for i=1:nMovies
    % set search path for parameter/intensity file
    if nargin>2 && ~isempty(channel)
        sourcePath = [getfield(data, {1,i}, channel) filesep];
    else
        sourcePath = [data(i).source filesep];
    end;
    
    if exist(sourcePath, 'dir')==7
        % get filenames for parameter files
        fileName_c0{i} = getParameterFileName(sourcePath, paramName_c0);
        fileName_c1{i} = getParameterFileName(sourcePath, paramName_c1);
        if (reference == 1)
            fileNameRef_c0{i} = getParameterFileName(sourcePath, paramNameRef_c0);
            fileNameRef_c1{i} = getParameterFileName(sourcePath, paramNameRef_c1);
        end;
        
        % if path is file name, use this file
    elseif exist(sourcePath, 'file')==2
        fileName_c0{i} = sourcePath;
    else
        error('No parameter file found.');
    end;
end;

%================================================================
% 2.
%================================================================
for i=1:nMovies
    fprintf('Movie no. %d\n', i);
    
    % extract the positions with the desired parameters
    % retrieves the indexes of the positions corresponding to the parameters (e.g. lifetime range from 'lftInfo.mat')
    %posvec = extractPos_vecFromLftRest([data(i).source filesep 'LifetimeInfo'], restvector, data(i).framerate);
    %nt = length(posvec); % number of CCPs in cohort
    %fprintf('Number of tracks in cohort: %d\n', nt);
    % load parameters stored as iMat_obj in parameterMat.mat
    paramMat_c0 = load(fileName_c0{i});
    paramMat_c0 = paramMat_c0.iMat_obj;
    paramMat_c0 = paramMat_c0(:,:,1); % extract intensities
    
    paramMat_c1 = load(fileName_c1{i});
    paramMat_c1 = paramMat_c1.iMat_obj;
    paramMat_c1 = paramMat_c1(:,:,1);
    
    nt = size(paramMat_c0, 1);
    
    % if reference parameter values are used, load background
    if (reference == 1)
        paramMatRef_c0 = load(fileNameRef_c0{i});
        paramMatRef_c0 = paramMatRef_c0.iMat_ref;
        
        paramMatRef_c1 = load(fileNameRef_c1{i});
        paramMatRef_c1 = paramMatRef_c1.iMat_ref;
        
        nbgPixels = size(paramMatRef_c0,1)/size(paramMat_c0, 1);
        % determine how many reference points exist for each data point
        % (e.g. 8 or 10 background points per data point)
        
        % if rfac>1 (e.g. rfac=8), then the assumption is that the 8 ref
        % points for data point number one are in the rows 1-8, those for
        % data point number 2 in rows 9-16, etc.
        % for easier use in the following function, the reference data are
        % immediately averaged here, i.e. they are compressed from rfac>1
        % to a matrix the same size as paramMat_use(:,:,1)
        
        bgMean_c0 = paramMatRef_c0(1:nbgPixels:end,:,1);
        bgMean_c1 = paramMatRef_c1(1:nbgPixels:end,:,1);
        for b = 2:nbgPixels
            bgMean_c0 = bgMean_c0 + paramMatRef_c0(b:nbgPixels:end,:,1);
            bgMean_c1 = bgMean_c1 + paramMatRef_c1(b:nbgPixels:end,:,1);
        end;
        bgMean_c0 = bgMean_c0 / nbgPixels;
        bgMean_c1 = bgMean_c1 / nbgPixels;
    end
    
    cargoStatus = zeros(1,nt);
    valid = zeros(1,nt);
    
    lambda = 20;
    for k = 1:nt
        %fprintf('k = %d\n', k);
        c0signal = paramMat_c0(k,:)';
        c0signal(isnan(c0signal)) = [];
        c1signal = paramMat_c1(k,:)';
        c1signal(isnan(c1signal)) = [];
        c0bg = bgMean_c0(k,:)';
        c0bg(isnan(c0bg)) = [];
        c1bg = bgMean_c1(k,:)';
        c1bg(isnan(c1bg)) = [];
        
        if ~isempty(c1bg)
            c1spline = smoothingSplineFourier(c1signal', lambda);
            c1bg_spline = smoothingSplineFourier(c1bg', lambda)';
            c1bg_std = sqrt(sum((c1bg-c1bg_spline).^2)/length(c1bg));
            
            omega = c0signal - c0bg;
            omega = omega / max(omega);
            
            w = (c1signal - c1bg).*omega + c1bg;
            w = w - c1bg_spline;
            pos = sum(w(w>0));
            neg = sum(w(w<0));
            
            %if (c1bg_spline(end) + 3*c1bg_std >= c1spline(end)) && (mean(c1bg) + std(c1bg) <= mean(c1signal))
            if (c1bg_spline(end) + 3*c1bg_std >= c1spline(end)) && (pos > -3*neg)
                cargoStatus(k) = 1;
            end;
            valid(k) = 1;
        end;
        %     figure;
        %     plot(c0signal, 'k');
        %     hold on;
        %     plot(c1signal, 'r');
        %     %plot(c0bg, 'k--');
        %     %plot(c1bg, 'r--');
        %     %plot(c1spline, 'k--');
        %     %plot(c1bg_spline, 'k');
        %     %plot(w, 'g');
        %     title([num2str(cargoStatus(k)) num2str(pos) num2str(neg)]);
        
    end;
    data(i).cargoStatus = cargoStatus;
    data(i).csValid = valid;
    save([data(i).source 'TrackInfoMatrices' filesep 'cargoStatus.mat'], 'cargoStatus', 'valid');
end;


function fileName = getParameterFileName(sourcePath, pfName)
if (exist([sourcePath  pfName], 'file')==2)
    fileName = [sourcePath pfName];
else
    [fileName, filePath] = uigetfile('*.mat', ['Select ' pfName], sourcePath);
    fileName = [filePath fileName];
end;


% function c = ncorr(a, b)
%
% a = a-mean(a(:));
% b = b-mean(b(:));
% c = sum(sum(sum(a.*b))) / sqrt(sum(sum(sum(a.^2)))*sum(sum(sum(b.^2))));


function [y y_dx y_d2x] = smoothingSplineFourier(s,lambda)

N = length(s);
s = [s s(N-1:-1:2)];
M = 2*N-2;

w = (0:M-1)*2*pi/M;

S = fft(s);

H = 3 ./ (2+cos(w)+6*lambda*(cos(2*w)-4*cos(w)+3));
c = real(ifft(S.*H));
c = c(1:N);

b3 = [1 4 1]/6;

y_dx = convn([c(2) c c(end-1)], [0.5 0.5 0], 'valid') - convn([c(2) c c(end-1)], [0 0.5 0.5], 'valid');
y_d2x = [c(2) c(1:end-1)] - 2*c + [c(2:end) c(end-1)];
y = convn([c(2) c c(end-1)], b3, 'valid');