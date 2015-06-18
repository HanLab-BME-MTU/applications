function [scaled_velocities] = scriptSTICS(imageStack)
%SCRIPTSTICS uses STICS analysis to output flow vectors between two images
% Short script to use STICS ImageJplug-in to determine flow between to images
% in a time lapse ideally
%Output:
% scaled_velocities is a 5xN array (N = number of vectors)
% 1st row: x velocity
% 2nd row: y velocity
% 3rd row: I honestly don't know, we don't use it
% 4th row: x position
% 5th row: y position

%Run Miji first to open ImageJ

% Make sure image is in 32-bit format and export to ImageJ
imageStack = uint32(imageStack);
MIJ.createImage(imageStack);

%Import necessary prefices
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jalgs.*;
import jalgs.jfft.*;
import jalgs.jseg.*;
import jguis.*;
import ij.measure.*;
import ij.text.*;

%Set parameters
imp = WindowManager.getCurrentImage();
MIJ.run('Properties...','slices=1 frames=10'); %Make sure image is read properly
MIJ.run('32-bit'); %Seriously, program ONLY runs this format
map = STICS_map(32,16);
stack=imp.getStack();
width=imp.getWidth();
height = imp.getHeight();
roi=Roi(0,0,width,height);

shift = 1;
frames = 10;
slices = 1;
channels = 1;
scaling = 1;
ftime = 1;
centered = true;
norm = true;
multiplier = 1;
magthresh = 2;
stepsize = 16;


%Update map with parameters and and run STICS analysis
map.update_STICS_map(jutils.get3DTSeries(stack,0,0,frames,slices,channels),width,height,0,frames,roi.getPolygon(),shift);
fp = map.get_map(scaling,ftime,stepsize,centered,norm,multiplier,stepsize,magthresh);
vector_stack=ImageStack(fp.getWidth(),fp.getHeight());
vector_stack.addSlice('',fp);
scaled_velocities=map.get_scaled_velocities(scaling,ftime,stepsize);
MIJ.run('Close All Without Saving')
end