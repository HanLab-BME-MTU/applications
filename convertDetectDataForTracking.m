function [movieInfo] = convertDetectDataForTracking(pathname, doubleStat)
% load detection data to be tracked
% INPUT: pathname =     location of the files to be read (assuming they are
%                       called cdata000x in a folder called maxdata283, which 
%                       is what Henry's software produces)
%       doubleStat =    parameter that determines what should happen to
%                       double maxima; for 
%                       =0 (default) the secondary maxima are read, provided
%                       they have a certain minimum distance from the
%                       primary maxima, and the intensities are split
%                       evenly
%                       =1 the secondary maxima are ignored, and primary
%                       maxima get the whole intensity  
% OUTPUT: movieInfo     = movieInfo is a file of the format that can be
%                       read and used for Khuloud's tracker
%
% Last modified: Francois Aguet, Feb 2010

if nargin < 2
    doubleStat = 0;
end

cdataList = dir([pathname filesep 'cdata*.mat']);
cdataList = {cdataList.name};
nFrames = length(cdataList);
movieInfo(1:nFrames) = struct('area', [], 'amp', [], 'xcoord', [], 'ycoord', []);

for k = 1:nFrames
    
    cdata = load([pathname filesep cdataList{k}]);
    % loading this file automatically loads:
    
    % csize = area
    % inn = 
    % intot = total intensity
    % lxm = number of maxima in object (1 or 2)
    % nmax = number of maxima
    % num = number of objects
    % xav = average x position (length num)
    % xmax = integer x-posititon of maximum (length nmax)
    % yav = average y position (length num)
    % ymax = integer y-posititon of maximum (length nmax)
    
    % if two maxima are found within one object, what happens?
    % possibility 1: the area and intensity are divided up evenly between
    % the individual objects, the missing (xav,yav) values are filled up
    % with the corresponding values of (xmax,ymax), which have integer
    % values
    % possibility 2: the secondary maximum is ignored
    % third possibility: accept only those secondary maxima that are
    % located more than a certain distance (dthresh) away from the primary
    % maximum
    dthresh = 5;
    
    useXcoord = cdata.xav; 
    useYcoord = cdata.yav;
    
    % for doubleStat=0 (default) : Version1 = splitting evenly
    if (doubleStat==0)
         
        ivec_posGoodMax = 1:cdata.num;
        ivec_numUseMax = ones(cdata.num,1);
        % if the number of maxima (nmax) is higher than the number of
        % discrete objects (num), the missing x,y values are filled up with
        % the xmax,ymax values; the secondary maxima are appended at the
        % end of the vector in cdata.xmax
        index_pm = zeros(cdata.num,1);
        interMaxDist = zeros(cdata.num,1);
        
        if (cdata.nmax>cdata.num)
            % loop over all secondary maxima:
            % find the positions of the corresponding primary maximum by
            % comparing intensities, and calculate the x,y-distance from
            % the secondary maximum
            for nm = cdata.num+1:cdata.nmax
                % position index of primary maximum corresponding to a
                % given secondary maximum
                fpos = find(cdata.intot == cdata.intot(nm));
                if length(fpos(fpos<=cdata.num))>2
                    display(['position of primary max corresponding to secondary max at pos ',num2str(nm),' is ambiguous']);
                end
                index_pm(nm) = min(fpos);
                % x,y av position of the primary maximum
                XYpos_pmAV = [cdata.xav(index_pm(nm)) cdata.yav(index_pm(nm))];
                % x,y max position of the primary maximum
                XYpos_pmMAX = [cdata.xmax(index_pm(nm)) cdata.ymax(index_pm(nm))];
                % x,y position of the secondary maximum
                XYpos_sm = [cdata.xmax(nm) cdata.ymax(nm)];
                % distance between primary and secondary maxima
                interMaxDist_av = sqrt((XYpos_pmAV(1)-XYpos_sm(1))^2+(XYpos_pmAV(2)-XYpos_sm(2))^2);
                interMaxDist_max = sqrt((XYpos_pmMAX(1)-XYpos_sm(1))^2+(XYpos_pmMAX(2)-XYpos_sm(2))^2);
                interMaxDist(nm) = min(interMaxDist_av,interMaxDist_max);
            end
            
            % if the distance interMaxDist is larger than the threshold
            % distance, include the secondary maximum in the list of
            % maxima, and also divide the intensity and area up evenly
            % between the primary and secondary maximum.
            goodPos = find(interMaxDist>dthresh);
            ivec_posGoodMax = [ivec_posGoodMax'; goodPos];
            
            if ~isempty(goodPos)
                % xy, positions of the usable secondary maxima are taken
                % from the xmax, ymax vectors
                useXcoord(cdata.num+1:cdata.num+length(goodPos)) = cdata.xmax(goodPos);
                useYcoord(cdata.num+1:cdata.num+length(goodPos)) = cdata.ymax(goodPos);
                % the index vector ivec_numUseMax is set to value+1 at the
                % position of the secondary AND of the corresponding
                % primary maximum, so that the measured intensity and area
                % can be divided evenly between the two (or three) maxima
                
                for smi = 1:length(goodPos)
                    % Primary Max positions - loop over these first to
                    % account for possible triple maxima
                    pmPos = index_pm(goodPos(smi));
                    ivec_numUseMax(pmPos)=ivec_numUseMax(pmPos)+1;
                end % for smi
                
                for smi = 1:length(goodPos)
                    % Secondary Max Positions
                    smPos = cdata.num+smi;
                    ivec_numUseMax(smPos)=ivec_numUseMax(index_pm(goodPos(smi)));
                end
            end % of if usable secondary maxima exist  
                      
        end

        % now the vector ivec_posGoodMax contains the index numbers of all
        % primary and usable secondary maxima, while the vector
        % ivec_numUseMax tells us whether the maximum at this position is a
        % single or one of a pair, and thus whether area and intensity are 
        % split evenly or not - we use this vector instead of cdata.lxm 
        % because it allows us to ignore those secondary maxima that don't
        % fit the criteria
        useArea = cdata.csize(ivec_posGoodMax)./ivec_numUseMax; 
        useInt = cdata.intot(ivec_posGoodMax)./ivec_numUseMax;        
    % else ignore the secondary maximum    
    else 
        % if the number of maxima (nmax) is higher than the number of
        % discrete objects (num), then the secondary maximum values are
        % simply ignored; area and intensity, which are duplicated in the 
        % vector at the position of the second maximum, are simply read out
        % only up to the number of objects
        useArea = cdata.csize(1:cdata.num);  
        useInt = cdata.intot(1:cdata.num); 
    end
    movieInfo(counter).area = useArea;
    movieInfo(counter).amp = useInt;
    movieInfo(counter).amp(:,2) = 0;
    
    movieInfo(counter).xCoord = useXcoord;
    movieInfo(counter).xCoord(:,2) = 0;
    
    movieInfo(counter).yCoord = useYcoord;
    movieInfo(counter).yCoord(:,2) = 0;
end