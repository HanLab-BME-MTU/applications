function estat = genSTORMtime(nEmi,kon,koff,kbleach,varargin)
%GENSTORMTIME generates time course of STORM dyes for calibration purposes
%   
%   Input:
%       nEmi     number of emitters
%       kon  on rate [1/sec] (from OFF to ON state)
%       koff     off rate [1/sec] (from ON to OFF)
%       kbleach  bleaching rate [1/sec]
%   optional:
%       nFrames  total number of frames, default 5000
%       expTime  exposure time [sec], default 0.1 sec
%
%   Output:
%       etime   time course for each emitter (0: off, 1: on, 2: bleached)
%
% US 11/07/2012
%

ip=inputParser;
ip.CaseSensitive=true;
ip.StructExpand=true;

ip.addRequired('nEmi',@isnumeric);
ip.addRequired('kon',@isnumeric);
ip.addRequired('koff',@isnumeric);
ip.addRequired('kbleach',@isnumeric);

ip.addParamValue('nFrames',5000);
ip.addParamValue('expTime',0.1);

ip.parse(nEmi,kon,koff,kbleach,varargin{:});

nEmi=ip.Results.nEmi;
kon=ip.Results.kon;
koff=ip.Results.koff;
kbleach=ip.Results.kbleach;
nFrames=ip.Results.nFrames;
expTime=ip.Results.expTime;

poff=koff/(koff+kbleach);

etime=zeros(nEmi,nFrames);

for e=1:nEmi
    xi=rand();
    if xi < kon*expTime
        etime(e,1)=1;
    end
end

for n=1:nFrames-1
    
    for e=1:nEmi

        switch etime(e,n)
            % emitter is in dark state
            case 0
                xi=rand();
                if xi < kon*expTime
                    etime(e,n+1)=1;
                end
                continue;
            % emitter is in fluorescent state
            case 1
                xi=rand();
                if xi < (koff+kbleach)*expTime
                    zeta=rand();
                    if zeta < poff
                        etime(e,n+1)=0;
                    else
                        etime(e,n+1)=2;
                    end
                else
                    etime(e,n+1)=1;
                end
                continue;
            % emitter is bleached    
            case 2
                etime(e,n+1)=2;
                continue;
        end
    end
end

% calculate number of blinking events, and ON, OFF and BLEACH time
estat=repmat(struct('etime',[],'nb',0,'ton',[],'toff',[],'tb',[]),nEmi,1);

for e=1:nEmi
    estat(e).etime=etime(e,:);

    if etime(e,1) == 1
        tonStart=1;
        tbStart=1;
    else
        tonStart=NaN;
        tbStart=NaN;
    end
    
    if etime(e,1) == 0
        toffStart=1;
    else
        toffStart=NaN;
    end
    
    tbEnd=NaN;
        
    for n=1:nFrames-1
        
        if ~isnan(tbEnd)
            break;
        end
        % emitter turns on
        if etime(e,n) == 0 && etime(e,n+1) == 1
            % increase number of blinks
            estat(e).nb=estat(e).nb+1;
            if isnan(tonStart)
                tonStart=n+1;
            end
            
            if isnan(tbStart)
                tbStart=n+1;
            end
            
            if ~isnan(toffStart)
                toffEnd=n;
                estat(e).toff=...
                    vertcat(estat(e).toff,[toffStart,toffEnd,toffEnd-toffStart+1]);
                toffStart=NaN;
            end
        end
        % emitter turns off
        if etime(e,n) == 1 && etime(e,n+1) == 0
            tonEnd=n;
            estat(e).ton=...
                vertcat(estat(e).ton,[tonStart, tonEnd, tonEnd-tonStart+1]);
            tonStart=NaN;
            
            if isnan(toffStart)
                toffStart=n+1;
            end
        end
        % emitter bleaches
        if etime(e,n) == 1 && etime(e,n+1) == 2
            if isnan(tbEnd)
                tbEnd=n;
            end
        end
        
    end
    
    if isnan(tbEnd)
        tbEnd=nFrames;
    end
    estat(e).tb=[tbStart, tbEnd, tbEnd-tbStart+1];
    estat(e).ton=vertcat(estat(e).ton,[tonStart,tbEnd,tbEnd-tonStart+1]);
            
end


end

