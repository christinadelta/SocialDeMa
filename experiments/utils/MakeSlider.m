function [set] = MakeSlider(scrn, set)

% sub-function to run the slider (scale) via beads task (confidence
% ratings), and phase 1 of best-choice tasks

windrect        = scrn.windrect;

% define variables 
anchors         = {'0', '50', '100'};
center          = round([windrect(3) windrect(4)]/2);
line            = 10;
width           = 3;
maxtime         = 10;       % maximum time for dragging the slider
scalelength     = 0.9;      % will change this?
scalepos        = 0.8;      % scale position (0 =top, 1=bottom, and in between)
startpos        = 'left';   % where will be the starting point of the slider?
mousebutton     = 1; 

% calculate coordinates of scale line 
midclick    = [center(1) windrect(4)*scalepos - line - 5 center(1), windrect(4)*scalepos + line + 5];
leftclick   = [windrect(3)*(1-scalelength) windrect(4)*scalepos - line windrect(3)*(1-scalelength) windrect(4)*scalepos + line];
rightclick  = [windrect(3)*(scalelength) windrect(4)*scalepos - line windrect(3)*(scalelength) windrect(4)*scalepos + line];
horzline    = [windrect(3)*scalelength windrect(4)*scalepos windrect(3)*(1-scalelength) windrect(4)*scalepos];

% Calculate the range of the scale, which will be need to calculate the
% position
scalerange          = round(windrect(3)*(1-scalelength)):round(windrect(3)*scalelength); % Calculates the range of the scale (0-100)
scalerangeshifted   = round((scalerange)-mean(scalerange)); % Shift the range of scale so it is symmetrical around zero

% UPDATE SETTINGS STRUCT 
set.anchors         = anchors;
set.scalelength     = scalelength;
set.mousebutton     = mousebutton;
set.leftclick       = leftclick;
set.rightclick      = rightclick;
set.midclick        = midclick;
set.horzline        = horzline;
set.scalerange      = scalerange;
set.width           = width;
set.scalepos        = scalepos;
set.startpos        = startpos;
set.line            = line;
set.maxtime         = maxtime;

end