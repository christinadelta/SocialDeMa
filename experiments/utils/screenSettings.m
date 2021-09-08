function [scrn] = screenSettings(scrn, taskNb)

% screen settings of the optimal stopping tasks
% 1. defines screen settings
% 2. defines object x and y 

%% define screen parameters 

% change these parameters as appropreate 
% screen resolution 
scrn.actscreenRes   = scrn.actscreen; % get screen's actual resolution
scrn.screenRes      = [1280 800]; % this also the windrect in px
scrn.hz             = 60; 
scrn.distview       = 700;
scrn.width          = scrn.actwidth;
scrn.height         = scrn.actheight;



%% define object x and y 
if taskNb == 3
    
    stimdeg             = scrn.stimdeg;
    % don't change anything here
    angleradn       = 2 * atan(scrn.width / 2 /scrn.distview);
    angledeg        = angleradn * (180/pi);
    pix             = scrn.screenRes(1) / angledeg;
    scrn.objectx    = round(stimdeg * pix);
    scrn.objecty    = round(stimdeg * pix);
    
end

end