function []=displayIntRes(intRes,align)

% this function displays the results gathered in the intRes structure
% INPUT:    intRes  = intensity results; this structure is expected to
%                   contain a number of fields (field names are
%                   irrelevant), and the intensity structures in subsequent
%                   fields are plotted
%           align (OPTIONAL): denotes whether alignment is to appearance or
%                   disappearance timepoint; if align==2, align to
%                   disappearance, otherwise, align to appearance

% color definitions

colors(1,:) = [0 0 1]; %blue
colors(2,:) = [0 1 1]; %cyan
colors(3,:) = [0 1 0]; %green
colors(4,:) = [1 0 0]; %red
colors(5,:) = [1 0 1]; %magenta;
colors(6,:) = [1 1 0]; %yellow
colors(7:12,:) = colors(1:6,:);
colors(13:18,:) = colors(1:6,:);

% loop over all fields in the structure (which can e.g. represent different lifetime cohorts)
 
cf = fieldnames(intRes); 

for n = 1:length(cf)

    % current field    
    cfield = intRes.(cf{n});
    
    % time vector
    if nargin>1 && align==2
        tvec = cfield.tvec(2,:);
    else
        tvec = cfield.tvec(1,:);
    end
   
    % intensity and error vectors
    ivec = nanmean(cfield.intAVE,1);
    evec = cfield.intSEM;

    % plot results
    errorbar(cfield.framerate*tvec,ivec,evec,'k-');
    hold on;
    plot( cfield.framerate*tvec,ivec,'Color', colors(n,:),'LineWidth',2 );
    
end

xlabel('time (sec)');
ylabel('intensity above BG');