%This script file shows the animation of identified body force over time.

%Load the identified body force.
load([resultPath 'bfId']);

M = vectorFieldAnimate([bfDisplayPx bfDisplayPy],recBF,10000, ...
   'bgImg',rawI{1},'vc','r');
