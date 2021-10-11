function [set] = MakeSlider(scrn, set)

% sub-function to make and draw the scale for rating via beads task (confidence
% ratings), and phase 1 of best-choice tasks

% UNPACK SCREEN RELATED STUFF
window          = scrn.window;      % main window
windrect        = scrn.windrect;
globalrect      = scrn.globalrect;
textsize        = scrn.textsize;
grey            = scrn.grey;
white           = scrn.white;
ycenter         = scrn.ycenter;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;

object_offset   = set.object_offset;
taskNb          = set.taskNb;
EEG             = set.EEG;

% define variables 
anchors         = {'0', '50', '100'};
center          = round([windrect(3) windrect(4)]/2);
line            = 10;
width           = 3;
scalelength     = 0.9;      % will change this?
scalepos        = 0.8;      % scale position (0 =top, 1=bottom, and in between)
startpos        = 'left';   % where will be the starting point of the slider?
mousebutton     = 1; 
scalecolour     = white;

if taskNb == 2
    thisprice   = set.thisprice;

elseif taskNb == 3
    textures    = set.textures;     % this will be used to draw the textures on screen
    thisitem    = set.thisitem;
    destrect    = set.destrect;
end

% UNPACK TRIGGER STUFF
if EEG == 1
    ioObj       = set.ioObject;
    status      = set.status;
    triggerdur  = set.triggerdur;
    address     = set.address;

    trigger9    = set.trigger11;
end

% First define the starting point of the slider
% calculate coordinates of scale line and text bounds
if strcmp(startpos, 'left')
   x = globalrect(3)*(1-scalelength);
elseif strcmp(startpos, 'center')
   x = globalrect(3)/2;
elseif strcmp(startpos, 'right')
   x = globalrect(3)*scalelength;
end

% ADD TEXTBOUNDS -- WILL BE USED TO CREATE THE SLIDER 
textbounds = [Screen('TextBounds', window, sprintf(anchors{1})); Screen('TextBounds', window, sprintf(anchors{3}))];    

% calculate coordinates of scale line 
midclick    = [center(1) windrect(4)*scalepos - line - 5 center(1), windrect(4)*scalepos + line + 5];
leftclick   = [windrect(3)*(1-scalelength) windrect(4)*scalepos - line windrect(3)*(1-scalelength) windrect(4)*scalepos + line];
rightclick  = [windrect(3)*(scalelength) windrect(4)*scalepos - line windrect(3)*(scalelength) windrect(4)*scalepos + line];
horzline    = [windrect(3)*scalelength windrect(4)*scalepos windrect(3)*(1-scalelength) windrect(4)*scalepos];

% Calculate the range of the scale, which will be need to calculate the
% position
scalerange          = round(windrect(3)*(1-scalelength)):round(windrect(3)*scalelength); % Calculates the range of the scale (0-100)
scalerangeshifted   = round((scalerange)-mean(scalerange)); % Shift the range of scale so it is symmetrical around zero

% CREATE WINDOWS FOR FLIPPING
if taskNb == 1
    rating_window = Screen('OpenOffscreenWindow',window);
    Screen('TextSize', rating_window, textsize);
    Screen('FillRect', rating_window, grey ,windrect);
    DrawFormattedText(rating_window, 'On a scale of 0 to 100, how confident are you for your choice?', 'center', ycenter-150, white);
    DrawFormattedText(rating_window, '0 = Not Confident', 'center', ycenter-50, white);
    DrawFormattedText(rating_window, '100 = Very Confident', 'center', ycenter, white);

elseif taskNb == 2
    
    rating_window = Screen('OpenOffscreenWindow',window);
    Screen('TextSize', rating_window, textsize);
    Screen('FillRect', rating_window, grey ,windrect);
    DrawFormattedText(rating_window, 'On a scale of 0 to 100, rate how likely it would be for you to get that contract in real life?', 'center', ycenter-200, white);
    DrawFormattedText(rating_window, '0 = I would never accept this contract', 'center', ycenter-150, white);
    DrawFormattedText(rating_window, '100 = I would definately accept this contract ', 'center', ycenter-100, white);

elseif taskNb == 3
    
    rating_window = Screen('OpenOffscreenWindow',window);
    Screen('TextSize', rating_window, textsize);
    Screen('FillRect', rating_window, grey ,windrect);
    DrawFormattedText(rating_window, 'On a scale of 0 to 100, rate how likely it would be for you date that person in real life.', 'center', ycenter-250, white);
    DrawFormattedText(rating_window, '0 = I would never date that person', 'center', ycenter-200, white);
    DrawFormattedText(rating_window, '100 = I would definately date that person', 'center', ycenter-150, white);
   
end

% Left, middle and right anchors
DrawFormattedText(rating_window, anchors{1}, leftclick(1, 1) - textbounds(1, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Left point
DrawFormattedText(rating_window, anchors{2}, 'center',  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Middle point
DrawFormattedText(rating_window, anchors{3}, rightclick(1, 1) - textbounds(2, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Right point

% Drawing the scale
Screen('DrawLine', rating_window, scalecolour, midclick(1), midclick(2), midclick(3), midclick(4), width);         % Mid tick
Screen('DrawLine', rating_window, scalecolour, leftclick(1), leftclick(2), leftclick(3), leftclick(4), width);     % Left tick
Screen('DrawLine', rating_window, scalecolour, rightclick(1), rightclick(2), rightclick(3), rightclick(4), width); % Right tick
Screen('DrawLine', rating_window, scalecolour, horzline(1), horzline(2), horzline(3), horzline(4), width);     % Horizontal line

% initialise the mouse
SetMouse(round(x), round(windrect(4)*scalepos), window, 1)

t0                         = GetSecs;
initresp                   = 0;

while initresp == 0 
    
    [x,~,buttons,~,~,~] = GetMouse(window, 1);
    
    % Stop at upper and lower bound
    if x > windrect(3)*scalelength
        x = windrect(3)*scalelength;
    elseif x < windrect(3)*(1-scalelength)
        x = windrect(3)*(1-scalelength);
    end
   
    Screen('CopyWindow', rating_window, window, windrect, windrect)
    
    if taskNb == 2
        DrawFormattedText(window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter, white); 
        
    elseif taskNb == 3
        Screen('DrawTexture', rating_window, textures{thisitem}, [], destrect); % display thisitem
        
    end
    
    object_onset = Screen('Flip', window, object_offset - slack);    % rating window is on
    
    % send confidence screen trigger
    if EEG == 1 
        io64(ioObj, address, trigger9)
        WaitSecs(triggerdur);
        io64(ioObj, address, 0) % return port to zero
    end
    
    % wait for second response 
    if buttons(mousebutton) == 1
        initresp = 1;
    end

end % end of while loop 

% release buttons
KbReleaseWait; %Keyboard

object_offset       = object_onset + 0.2 - ifi;
set.object_offset   = object_offset;
set.object_onset    = object_onset;

return