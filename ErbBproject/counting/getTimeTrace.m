function [traces,fullTraces,estat,erg]=getTimeTrace(pList,stack,varargin)
%GETTIMETRACE calculates fluorescent time traces of single molecules
%
%   required input arguments:
%       pList -> particle list in each bin
%       stack -> image stack
%
%   optional input arguments:
%          nRep -> number of frames within activation cycle, default 10
%       display -> display results: true/false, default 'false'
%
%   output:
%           traces -> cell array with time traces
%       fullTraces -> matrix with full time trace, from beginning to end
%           estate -> fluorescent state of detected molecule, ON 1, OFF 0
%              erg -> structure, exponential fit of time scales to:
%                     time spent in ON state  ->  .onScale
%                     time spent n OFF state  ->  .offScale
%                     time until bleaching    ->  .bScale
%                     number of activations   ->  .actScale
%
%   US 2012/11/21
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('pList',@iscell);
ip.addRequired('stack',@isnumeric);
ip.addOptional('nRep',10,@isscalar);
ip.addOptional('display',false,@islogical);

ip.parse(pList,stack,varargin{:});

pList=ip.Results.pList;
stack=ip.Results.stack;
nRep=ip.Results.nRep;

nFrames=size(stack,3);

np=numel(pList);
traces=cell(np,1);
fullTraces=zeros(np,nFrames);

for k=1:np
    info=pList{k};
    xx=ceil(info(1,2));
    yy=ceil(info(1,1));
    last=size(info,2);
    ff=info(:,last);
    traces{k}=squeeze(stack(xx,yy,:));
    
    fullTraces(k,ff)=1;
end

% calculate number of blinking events, and ON, OFF and BLEACH time
estat=repmat(struct('nAct',0,'tAct',[],'ton',[],'toff',[],'tb',[]),np,1);

% frames directly after UV activation:
idAct=1:nRep:nFrames;

for k=1:np

    % activation times and frequency of activation of each emitter
    allAct=fullTraces(k,idAct);
    nAct=sum(allAct);
    tAct=allAct.*idAct;
    tAct=tAct(tAct > 0);
    
    estat(k).nAct=nAct;
    estat(k).tAct=tAct;
    
    % length of active state in frames
    ton=NaN(nAct,1);
    for i=1:nAct
        t=tAct(i);
        on=fullTraces(k,t);
        count=0;
        while on
            count=count+1;
            if t+count <= nFrames
                on=fullTraces(k,t+count);
            else
                on=0;
            end
        end
        ton(i)=count;
    end
    estat(k).ton=ton;
    
    % length of inactive state in frames
    toff=NaN(nAct-1,1);
    for i=1:nAct-1
        toff(i)=tAct(i+1)-tAct(i);
    end
    estat(k).toff=toff;
    
    % time until photbleaching occurs
    %tb=tAct(end)-tAct(1);
    if ~isempty(tAct)
        estat(k).tb=tAct(end);
    end
end

% estimate scale of exponential distribution
ton=vertcat(estat.ton);
[muhat,muci]=expfit(ton);
erg.onScale=[muhat muci'];

toff=vertcat(estat.toff);
[muhat,muci]=expfit(toff);
erg.offScale=[muhat muci'];

tb=vertcat(estat.tb);
[muhat,muci]=expfit(tb);
erg.bScale=[muhat muci'];

nAct=vertcat(estat.nAct);
[muhat,muci]=expfit(nAct);
erg.actScale=[muhat muci'];

if ip.Results.display
    % display results
    tcolor=mat2gray(squeeze(label2rgb(1:np)),[0,255]);
    
    hold on
    for k = 1:np
        plot(traces{k},'-','Color',tcolor(k,:));
    end
end

end
    
    

