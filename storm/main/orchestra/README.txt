HOW TO RUN JOBS ON ORCHESTRA
============================

1) Check out the "qfsm" and the "common" repository on orchestra

2) Make sure the MATLAB path is set up correctly and all the workers 
   can access the files in the repositories.

3) Choose a location on a shared drive where you want to store all the 
   STORM related data and create the following folder structure.
   
   STORM-ROOT
        |-----------"_data"
        |-----------"_queue"
        |-----------"_presets"
        |-----------"_deleted"
        |-----------"_archive"

4) Set the "STORM-ROOT" data path in qfsm/storm/main/stormPath.m

5) Place your data sets in "_data" in a subfolder. Ideally it contains
   a "*.bin" and a "*.ini" file.

6) Use guiMain.m to create a job configuration and label them with 
   the "-" to run them in batch mode on orchestra.

8) Run one of the provided scripts (E.g. silentMax.sh) to process the 
   labeled data sets. If you use the "silentMax.sh" script you should 
   receive an email once the jobs are finished. 