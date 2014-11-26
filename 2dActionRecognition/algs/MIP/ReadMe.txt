
----------------------- MIP Implementation --------------------

Last updated:  25/10/12
 
http://www.openu.ac.il/home/hassner/projects/MIP/MIPcode.zip

=======
Authors
=======
 
The MIP (Motion Interchange Patterns) descriptor was developed in 2012
by Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf,
as part of an Action Recognition system.
For further details on the system, applications and experiments please refer to the following paper:

Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
Motion Interchange Patterns for Action Recognition in Unconstrained Videos,
European Conference on Computer Vision (ECCV), Firenze, Italy, Oct 2012.


The SOFTWARE ("MIPcode") is provided "as is", without any guarantee made
as to its suitability or fitness for any particular use. It may contain bugs, so use
of this tool is at your own risk. We take no responsibility for any damage that
may unintentionally be caused through its use.
 
For questions and bug reports please contact Tal Hassner (hassner@openu.ac.il)
 
=======
General
=======
 
The code in this directory computes the Motion Interchange Patterns descriptor
for a given video or all videos in set of videos.
The code is a slightly different version than the one described in the ECCV'12 paper
in the sense that it preforms alignment once for all channels.
This leads to better run-time performance and appears to produce similar results.
 
In this ReadMe file we specify the MIP code usage through a very simple GUI.
You can run this GUI both on Unix and on Windows (parfor is supported).

Otherwise, if you wish to run this code directly from Matlab's command line:
The main function is: run_MIP_coding_on_files (test_dir,params,database_info,add_path)
Please see the requested input at the function header.

Dependencies:
=============
To run this code you should have Matlab 2010 or higher,
This is required for Matlab's "mmreader" video input function.
You should also have the appropriate video codec installed on your system.

For debug mode - you should also install the following codec on your computer
http://www.xvid.org/Downloads.15.0.html
(NOTE: the default debug option is 'false', so no debug files will be written.)

=====
Files
=====
 
The MIPcode directory is divided into two sub-directories:
 
Code & Data
===========
The 'Code' directory includes the following:
- GUI               - The GUI code
- Features      - The actual Matlab code
 
The 'Data' directory includes the following:
- InputDir      - default directory for the videos you would like to code.
- OutputDir   - default directory for the current test output descriptors.
These directories currently include examples from the ASLAN set.
 
=======
Running
=======
 
Basic GUI instructions
----------------------------
 
1. Run MATLAB.
 
2. Go to the GUI directory (MIPcode/Code/GUI)
The GUI must run from this directory!
 
3. Enter: MIP_Coding_GUI on the command line.
 
4. Define the Movies Input Dir by browsing (default is MIPcode/Data/InputDir).
This directory should point to where the videos you would like to code.
The program will search for all the video files with a given extension (default: .avi and .wmv)
under this directory and all its sub-directories,
and will give the list of these files under 'Files List'.
From these files you can select to code one, a few, or all files by
clicking shift/ctrl + mouse left click
then choose either: All Files (default) / Selected Files
on the right side of the GUI.
 
5. Define an Out Dir (default is  MIPcode/Data/OutputDir).
This is the current test directory and where the output descriptors will be written.
 
6. To change parameters - click on:
'Change Params' on the upper left side of the GUI.
This will create and open the file:
'get_feature_params_runtime.m' as a text file,
and you can change the parameters there.
NOTE: This will only change the parameters of the current run,
and will NOT change the default parameters of other runs.

For example: to change the debug option,
you should have: params.MIP_params.MIP_dbg_dumps = true;
and save the 'fixed get_feature_params_runtime.m' file.
 
=======================================================================
If you use any of our data or code we ask that you cite the following:
=======================================================================
 
Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf,
Motion Interchange Patterns for Action Recognition in Unconstrained Videos,
European Conference on Computer Vision (ECCV), Firenze, Italy, Oct 2012.
 
BibTex entry:
-----------------
@inproceedings{K-GGHW:ECCV12:MIP,
 author    = {Orit~Kliper-Gross and Yaron~Gurovich and Tal~Hassner and Lior~Wolf},
 title     = {Motion Interchange Patterns for Action Recognition in Unconstrained Videos},
 booktitle = {European Conference on Computer Vision ({ECCV})},
 month  =  {Oct.},
 year   = {2012},
 URL    = {http://www.openu.ac.il/home/hassner/projects/MIP}
}
 
Questions and comments can be sent to:
----------------------------------------------------
orit.kliper@weizmann.ac.il
yarongur@post.tau.ac.il
hasnner@openU.ac.il
wolf@tau.ac.il
 
Thank you,
The authors.


