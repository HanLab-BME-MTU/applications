function movieData = calculateProtrusionAutocorrelation(movieData)


p.BatchMode = false;
p.MaxLag = 200;


%UNDER CONSTRUCTION

iProtProc = movieData.getProcessIndex('ProtrusionProcess',1,~p.BatchMode);

if ~movieData.processes_{iProtProc}.checkChannelOutput
    error('Movie must have valid protrusion vectors calculated first!')
end

tmp = movieData.processes_{iProtProc}.loadChannelOutput;

protrusion = tmp.protrusion;
smoothedEdge = tmp.smoothedEdge;
normals = tmp.normals;

nFrames = movieData.nFrames_;

nPtsTot = 0;

for iFrame = 1:(nFrames-1)
    

    %Normalize the normals - sam's function does not return them with unit
    %length
    normals{iFrame} = normals{iFrame} ./ ...
                repmat(sqrt(dot(normals{iFrame}',normals{iFrame}')'),1,2);
    

    %Get the normal component of the protrusion at each point
    ncCurr = dot(normals{iFrame}',protrusion{iFrame}')';
    %Get the magnitude of the protrusion at each point
    mCurr = sqrt(dot(protrusion{iFrame}',protrusion{iFrame}'))';

    %use Khulouds splitting trick - divide more than once??? Why only two??
    %Middle justified? Or something fancier?    
    nPtsCurr = numel(mCurr);
    nHalf = round(nPtsCurr/2);    
    normalComponent(iFrame).observations = ncCurr(1:nHalf);
    normalComponent(iFrame+nFrames-1).observations = ncCurr(nHalf+1:end);        
    magnitude(iFrame).observations = mCurr(1:nHalf);
    magnitude(iFrame+nFrames-1).observations = mCurr(nHalf+1:end);
    nPtsTot = nPtsTot + nPtsCurr;
    
end


%Get the total autocorrelation
[acNormalComponent,err(1)] = autoCorr(normalComponent,p.MaxLag);
[acMagnitude,err(2)] = autoCorr(magnitude,p.MaxLag);
    
if p.BatchMode
    acFig = figure('Visible','off');
else
    acFig = figure;
end

plot(0:p.MaxLag,acNormalComponent(:,1),'r');
hold on
plot(0:p.MaxLag,acMagnitude(:,1),'b');
plot(xlim,ones(1,2)*1.96/sqrt(nPtsTot),'--k');
legend('Normal Component','Magnitude','+/-1.96/sqrt(n)');
plot(xlim,-ones(1,2)*1.96/sqrt(nPtsTot),'--k');
plotTransparent(0:p.MaxLag,acNormalComponent(:,1),acNormalComponent(:,2),'r',.2,0);
plotTransparent(0:p.MaxLag,acMagnitude(:,1),acMagnitude(:,2),'b',.2,0);
if any(err)
    error('Problem calculating autocorrelation!')
end
