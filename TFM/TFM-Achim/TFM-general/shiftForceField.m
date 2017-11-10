function [forceField]=shiftForceField(forceField,displField,frameList)
% This function performs a shift of the force field given the underlying
% displacement field. Two analysis are performed:
%
% 1) shift of only the positions pos by the displ. at these pos.
%
% 2) reinterpolation of the shifted force field on the resulting
%    non-regular grid on the initial (ususally regular) grid positions.
%
% While 1) can be useful in all projects, the field resulting from 2) seems
% to provide no real improvement. The forces at the cell boundary seem to
% be persistently smaller, the forces within the cell are pretty much the
% same. Thus it might be that the reinterpolation step causes a systematic
% underestimation of the forces.

if nargin < 1 || isempty(forceField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select force field to be used as overlay');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       forceField=fileStruct.forceField;
end

if nargin < 2 || isempty(displField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select displField to be used for shifting the forceField');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       displField=fileStruct.displField;
end

%get the target directory:
if nargin < 3 || isempty(frameList)
    frameList=1:min(length(forceField),length(displField));
end

for i=frameList
    display(['Force Field of frame: ',num2str(i),' is shift according to the measured displacement field and then re-interpolated onto the initial grid positions.']);
    % interpolate the displacement field on to the grid of the forceField.
    ux_intp = TriScatteredInterp(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1));
    uy_intp = TriScatteredInterp(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,2));
    
    ux_at_fpos = ux_intp(forceField(i).pos(:,1),forceField(i).pos(:,2));
    uy_at_fpos = uy_intp(forceField(i).pos(:,1),forceField(i).pos(:,2));
    
    % some times there is an extra-polation at the grid boundary which
    % generates NaN and causes problems later on. We set those values to 0.
    ux_at_fpos(isnan(ux_at_fpos))=0;
    uy_at_fpos(isnan(uy_at_fpos))=0;
    
    % now add this to the position of the force vectors:
    posShifted(:,1)=forceField(i).pos(:,1)+ux_at_fpos;
    posShifted(:,2)=forceField(i).pos(:,2)+uy_at_fpos;
    
    % store this values in the force field structure:
    forceField(i).posShifted=posShifted;
    
    % these remain unchanged:
    vecShifted=forceField(i).vec;
    
    % Now reinterpolate the force field onto the initial grid positions:
    fx_shifted = TriScatteredInterp(posShifted(:,1),posShifted(:,2),vecShifted(:,1));
    fy_shifted = TriScatteredInterp(posShifted(:,1),posShifted(:,2),vecShifted(:,2));
    
    fx_ReIntp = fx_shifted(forceField(i).pos(:,1),forceField(i).pos(:,2));
    fy_ReIntp = fy_shifted(forceField(i).pos(:,1),forceField(i).pos(:,2));
       
    % Sometimes there is an extra-polation at the grid boundary which
    % generates NaN. Here we have to keep these values. Set the new values:
    forceField(i).vecReIntp(:,1)=fx_ReIntp;
    forceField(i).vecReIntp(:,2)=fy_ReIntp;
    
    clear posShifted posShifted
end 
    
    
    
    
    