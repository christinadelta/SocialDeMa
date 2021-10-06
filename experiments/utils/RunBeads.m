function [set, logs] = RunBeads(set, scrn, logs)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
taskNb          = set.taskNb;       % number of task (needed to run the correct task)
thistrial       = set.thistrial;    % this is the current sequence/trial number 
thisblock       = set.thisblock;     % this is the current block number 

isi             = set.isi;          % interstimulus interval
jitter          = set.jitter; 
fixation        = set.fixation;     % draw fixation
EEG             = set.EEG;          % is it a EEG session?
triggerdur      = set.triggerdur;   % trigger duration (3 ms)

window          = scrn.window;       % main window
windrect        = scrn.windrect;
xcenter         = scrn.xcenter;
ycenter         = scrn.ycenter;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;        % slack is ifi/2 (important for timing)
white           = scrn.white;
grey            = scrn.grey;
fixsize         = scrn.fixationsize;
textfont        = scrn.textfont;
textsize        = scrn.textsize;
textbold        = scrn.textbold;
smalltext       = scrn.smalltext;
baserect        = scrn.baserect;
centeredrect    = scrn.centeredrect;
penwidth        = scrn.penwidth;

thisession      = logs.sess;
sub             = logs.sub;
taskname        = logs.task;
resfolder       = logs.resultsfolder;

xcntr           = xcenter - 30; % works better than 'center'

% create fixation cross offscreen and paste later (faster)
fixationdisplay = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixationdisplay, grey);
Screen('TextFont',fixationdisplay, textfont);
Screen('TextSize',fixationdisplay, fixsize);
DrawFormattedText(fixationdisplay, fixation, 'center', ycenter, white);

%% ----- Run the beads task ------ %%

abort           = 0;
HideCursor;

% UNPACK BEADS RELATED PARAMETERS 
green           = scrn.green;
blue            = scrn.blue;
red             = scrn.red;
orange          = scrn.orange;

sequence        = set.sequence;     % this is the current sequence/trial 
thisurn         = set.urn;          % this is the current urn, contains the high prob colour beads 1=blue, 0=green
drawlen         = set.draws;        % length of each sequence (up to 10 draws)
penalty         = set.penalty;      % £0.25 loss for every draw
win             = set.win;          % £10 for every win!
balance         = set.balance;      % balance starts as zero and updates from correct/incorrect (loss) responses & draw choices

bead_dur        = set.bead_dur;     % duration of each bead in sec
response        = set.response;     % response duration in sec
fix_dur         = set.fix_dur;      % fixation duration
dfeedback       = set.feed_dur;     % duration of feedback window in sec

bluekey         = set.code1;       % keycode for blue option
greenkey        = set.code2;       % keycode for green option
drawkey         = set.code3;       % keycode for draw-again option
esckey          = set.code21;      % keycode for aborting the experiment 

% ADD A FEW MORE VARIABLES
blocktrials     = []; % store trial info
draws           = []; % store info specifically about the sequence/draws [1-10]
previous        = []; % store the previous draws here to show at the bottom of the screen
colours         = []; % store the previous draw-colours here to show at the bottom of the screen
draw_count      = 0;  % drawing counter
accuracy        = nan;

% UNPACK EEG TRIGGERS
if EEG == 1 
    
    ioObj       = set.ioObject;
    status      = set.status;
    
    trigger1    = set.trigger1;
    trigger2    = set.trigger2;
    trigger3    = set.trigger3;
    trigger4    = set.trigger4;
    trigger5    = set.trigger5;
    
    trigger6    = set.trigger6;
    trigger7    = set.trigger7;
    trigger8    = set.trigger8;
    trigger14   = set.trigger14;
    trigger15   = set.trigger15;
    trigger16   = set.trigger16;

    
    trigger102  = set.trigger102; % start of sequence
    trigger103  = set.trigger103; % end of sequence
    
end
    
% SET ALL RELEVANT WINDOWS (BEAD WINDOWS, RESPONSE WINDOW, CONFIDENCE RATING & FEEDBACK WINDOWS)
% BASE RECTANGLES -- white rect
whiterect_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', whiterect_window, textsize);
Screen('FillRect', whiterect_window, grey ,windrect);
Screen('FrameRect', whiterect_window, white, centeredrect, penwidth)

% BASE RECTANGLES -- orange rect
orangerect_window = Screen('OpenOffscreenWindow',window);
Screen('FillRect', orangerect_window, grey ,windrect);
Screen('FrameRect', orangerect_window, orange, centeredrect, penwidth);

% FEEDBACK WINDOW (you win!)
feedback_window1 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window1, textsize);
Screen('FillRect', feedback_window1, grey ,windrect);
DrawFormattedText(feedback_window1, 'Well Done! :)', 'center', ycenter, green);

% FEEDBACK WINDOW (you lose!)
feedback_window2 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window2, textsize);
Screen('FillRect', feedback_window2, grey ,windrect);
DrawFormattedText(feedback_window2, 'Oh No! :(', 'center', ycenter, red);

% FEEDBACK WINDOW (wrong answer!)
feedback_window3 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window3, textsize);
Screen('FillRect', feedback_window3, grey ,windrect);
DrawFormattedText(feedback_window3, 'Oh No! :(', xcntr, ycenter-50, red);
DrawFormattedText(feedback_window3, 'You you are not allowed to draw again', 'center', ycenter, red);

%%
% START THE TRIAL WITH A FIXATION CROSS
Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
fliptime    = Screen('Flip', window); % flip fixation window
trialstart  = fliptime;

% send sequence start trigger
if EEG == 1 
    io64(ioObj, address, trigger102)
    WaitSecs(triggerdur);
    io64(ioObj, address, 0)
end

% object offset
object_offset   = trialstart + fix_dur + isi + randperm(jitter*1000,1)/1000 - ifi;

% 1. BEGIN DRAWING 
for thisdraw = 1:drawlen
    
    % allow subject to abort the experiment 
    [keyisdown,~,keycode] = KbCheck;
    if keyisdown && keycode(esckey)
        abort = 1;
        break;
    end
    
    xcenters        = xcntr - 70; % start adding the darws prices at xcenter - 30
    xs              = xcenters;
    
    % 1. SHOW RECT WITH THE LIST OF BEADS(S)
    Screen('CopyWindow', whiterect_window,window, windrect, windrect)
    Screen('TextSize', window, textsize);
    Screen('TextStyle', window, textbold)

    % COND 1: Is it a high or a low probability draw?
    if sequence(thisdraw) == 1  % high prob draw

        % COND 2: Is it a blue or a green bead?
        if thisurn == 1         % blue urn and blue bead

            % show blue window
            DrawFormattedText(window, 'BLUE', 'center', ycenter, blue);
            if EEG == 1
                stimtrigger         = trigger1; % assign trigger 1 as stimulus triger
            end

            previous_colour     = blue;
            previous_bead       = 'BLUE';

        else % green urn and green beed

            % show green window
            DrawFormattedText(window, 'GREEN', 'center', ycenter, green);
            if EEG == 1
                stimtrigger         = trigger3;
            end

            previous_colour     = green;
            previous_bead       = 'GREEN';

        end % end of cond 2 loop

    elseif sequence(thisdraw) == 2 % low prob draw

        if thisurn == 1 % blue urn and green bead

            DrawFormattedText(window, 'GREEN', 'center', ycenter, green);
            if EEG == 1
                stimtrigger         = trigger2;
            end
            previous_colour     = green;
            previous_bead       = 'GREEN';

        else % % green urn and blue bead

            DrawFormattedText(window, 'BLUE', 'center', ycenter, blue);
            if EEG == 1
                stimtrigger         = trigger4;
            end

            previous_colour     = blue;
            previous_bead       = 'BLUE';

        end % % end of cond 2 loop
    end % end of cond 1 loop
    
    if thisdraw > 1
        
        previous_len        = length(previous); % how many previous prices?
        
        for i = 1:previous_len
            
            % show the previous prices on the left of the current bead
            Screen('TextSize', window, smalltext);
            Screen('TextStyle', window, 0)
            DrawFormattedText(window, previous{i}, xcenters, ycenter, colours{i}); 

            % update xcenters, so that previous prices are not
            % displayed on top of each other
            xcenters = xcenters - 60;
            
        end  
    end
    bead_onset     = Screen('Flip', window, object_offset - slack); 
    bead_offset    = bead_onset + bead_dur - ifi;                          % bead on for 0.5 ms
    
    % send sequence start trigger
%     if EEG == 1
%         sp.sendTrigger(stimtrigger) % blue urn -- high prob blue bead trigger
%     end
    
    % 2. SHOW GO SIGNAL
    % SHOW RECT WITH THE LIST OF BEADS(S)
    Screen('CopyWindow', orangerect_window,window, windrect, windrect)
    Screen('TextSize', window, smalltext);
    Screen('TextStyle', window, 0) % NOT BOLD
    DrawFormattedText(window, previous_bead, 'center', ycenter, previous_colour);
    Screen('TextSize', window, textsize);
    DrawFormattedText(window, '1:', xcntr-180, ycenter+220, orange);
    DrawFormattedText(window, 'blue?', xcntr-200, ycenter+250, orange);
    DrawFormattedText(window, '2:', 'center', ycenter+220, orange);
    DrawFormattedText(window, 'green?', 'center', ycenter+250, orange);
    DrawFormattedText(window, '3:', xcntr+200, ycenter+220, orange);
    DrawFormattedText(window, 'draw again?', xcntr+160, ycenter+250, orange);
    
    if thisdraw == drawlen
        Screen('TextSize', window, textsize);
        DrawFormattedText(window, 'This was your last draw.', 'center', ycenter-200, white);
    end  
    
    if thisdraw > 1 
        
        previous_len = length(previous);
        
        for i = 1:previous_len
            
            % show the previous prices on the left of the current bead
            Screen('TextSize', window, smalltext);
            Screen('TextStyle', window, 0)
            DrawFormattedText(window, previous{i}, xs, ycenter, colours{i}); 

            % update xcenters, so that previous prices are not
            % displayed on top of each other
            xs = xs - 60;
        end
    end
        
    response_onset    = Screen('Flip', window, bead_offset - slack); 
    
%     % send response prompt trigger
%     if EEG == 1 
%         sp.sendTrigger(trigger5)
%     end
    
    % 3. INIT RESPONSE: initiate on each draw
    rt                  = NaN;
    answer              = NaN;
    resp_input          = 0;
    responseTrigNotSent = 1;
    
    while resp_input == 0 && (GetSecs - response_onset) < response 
        
        [keyisdown, secs, keycode] = KbCheck;
        pressedKeys = find(keycode);
        
        % check the response key
        if isempty(pressedKeys) % if something weird happens or subject doesn't press any of the correct keys
            resp_input  = 0; 
            rt          = nan;
            answer      = nan;
            respmade    = GetSecs;
            
        elseif ~isempty(pressedKeys) % if subject pressed a valid key
            
            if keycode(1,bluekey) % if subject chose the blue urn 
                resp_input  = bluekey;
                rt          = secs - response_onset;
                answer      = 1; % blue urn 
                respmade    = secs;
                
%                 % send blue choice trigger
%                 if EEG == 1 && responseTrigNotSent==1
%                     sp.sendTrigger(trigger6);
%                     responseTrigNotSent=0;
%                 end
                
            elseif keycode(1,greenkey) % if subject chose the green urn 
                resp_input  = greenkey;
                rt          = secs - response_onset;
                answer      = 2; % green urn
                respmade    = secs;
                
%                 % send green choice trigger
%                 if EEG == 1 && responseTrigNotSent==1
%                     sp.sendTrigger(trigger7);
%                     responseTrigNotSent=0;
%                 end
                
            elseif keycode(1,drawkey) % if subject chose to draw again
                resp_input  = drawkey;
                rt          = secs - response_onset;
                answer      = 3; % draw-again 
                respmade    = secs;
                draw_count  = draw_count + 1;
                
%                 % send draw again trigger
%                 if EEG == 1 && responseTrigNotSent==1
%                     sp.sendTrigger(trigger8);
%                     responseTrigNotSent=0;
%                 end
            end %  
        end % if responded statement
    end % end of response while loop 
    
    prompt_offset   = respmade + isi - ifi;                                     % response prompt self paced or on for 2500 ms
    
   % 4. BRING BACK FIXATION. NEXT STEP:
   % A) if answer = 3, draw again,
   % B) if answer = 1 or answer = 2, ask subject to rate the confidence of their choice,
   % check if response was correct or incorrect, give appropreate feedback, break 
   % C) if draw_count == 10 and subject chose to draw again collect it as
   % wrong answer 
   Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
   fixation_onset   = Screen('Flip', window, prompt_offset - slack);                % fixation on, prepare for next draw or for feedback      
   object_offset    = fixation_onset + fix_dur + isi + randperm(jitter*1000,1)/1000 - ifi;     % add jitter here?
   
   if answer ~= 3
       
       if any(isnan(answer)) % if subject doesn't respond, just move to the next draw
           
           % if this is the last darw and subject didn't respond, update
           % accuracy and show "you lose" screen
           if thisdraw == drawlen
               accuracy = 0;
               
               Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
               feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 

               object_offset = feedback_onset + dfeedback - ifi;

               % send you lose trigger
%                if EEG == 1
%                    sp.sendTrigger(trigger15);
%                end
               
               % BRING FIXATION BACK ON
               Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
               object_onset     = Screen('Flip', window, object_offset - slack); 

               object_offset    = object_onset + isi - ifi;% no jitter
               
           else
               draw_count  = draw_count + 1; % if the subject didn't respond and this is not the last darw, judt update draw count and show next draw
           end
    
       else
           
           set.object_offset = object_offset; % 
           
           % THIS IS THE PART WHERE WE DRAW THE SLIDER AND SHOW THE
           % CONFIDENCE SCREEN
           set = MakeSlider(scrn, set); % first draw the scale and allow subject to click the mouse
           WaitSecs(0.1)                % delay to prevent fast mouse clicks mix 

           set = RunSlider(scrn, set);  % once subject click the first time, display slider
           WaitSecs(0.1)
               
           % UNPACK SETTINGS
           object_offset    = set.object_offset;
           rate_rt          = set.rate_rt;
           position         = set.position;

           if answer == 1 % if subject chose blue urn
               if thisurn == 1                                                      % if thisurn is in fact a blue urn

                   accuracy         = 1;                                            % update the accuracy var

                   Screen('CopyWindow', feedback_window1,window, windrect, windrect)% show "you win" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);  % feedback window on 

                   % send you win trigger
%                     if EEG == 1
%                         sp.sendTrigger(trigger14);
%                     end

               else % if "thisurn" is a green urn 
                   accuracy         = 0;

                   Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 

                   % send you lose trigger
%                     if EEG == 1
%                         sp.sendTrigger(trigger15);
%                     end
               end
               object_offset = feedback_onset + dfeedback - ifi;
               
               % show fixation between ending the trial/sequence
               Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
               object_onset     = Screen('Flip', window, object_offset - slack); 

               object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
               break;

           elseif answer == 2 % if subject chose green urn 

               if thisurn == 1 % if thisurn is a blue urn

                   accuracy         = 0;

                   Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 
                   % send you lose trigger
%                     if EEG == 1
%                         sp.sendTrigger(trigger15);
%                     end

               else % if "thisurn" is a green urn 

                   Screen('CopyWindow', feedback_window1,window, windrect, windrect) % show "you win" feedback window
                   accuracy = 1;                                                     % update the accuracy var
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 
                   % send you lose trigger
%                    if EEG == 1
%                        sp.sendTrigger(trigger14);
%                    end
               end
               object_offset = feedback_onset + dfeedback - ifi;
               
               % show fixation between ending the trial/sequence
               Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
               object_onset     = Screen('Flip', window, object_offset - slack); 

               object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
               break

           end
       end
       
   elseif answer == 3
       
       if thisdraw == drawlen % if this is the last (10th) draw ans subject requested another draw
           
           accuracy     = 0;
           % DISPLAY THE 3RD FEEDBACK SCREEN
           Screen('CopyWindow', feedback_window3,window, windrect, windrect)
           feedback_onset   = Screen('Flip', window, object_offset - slack); % feedback window on 
           % send you lose (out of draws) trigger
%            if EEG == 1
%                sp.sendTrigger(trigger16);
%            end
           
           % update answer, given that participant pressed 3 (for draw 10), answer should be 0 (this is an incorrect trial)
           answer           = 0;
           
           object_offset    = feedback_onset + dfeedback - ifi;
           
           % BRING FIXATION BACK ON BEFORE ENDING THE TRIAL
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           object_onset     = Screen('Flip', window, object_offset - slack); 
           
           object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
           break
       end
   end
   
   % store the current price to show at the bottom of the screen
   % (during the next samples)
   previous{thisdraw}          = previous_bead;
   colours{thisdraw}           = previous_colour;
           
   % add sequence/draws related info here 
   draws(thisdraw).session      = thisession;
   draws(thisdraw).block        = thisblock;
   draws(thisdraw).trialnumber  = thistrial;
   draws(thisdraw).trialonset   = trialstart;
   draws(thisdraw).beadonset    = bead_onset;
   draws(thisdraw).thisdraw     = thisdraw;
   draws(thisdraw).thisbead     = sequence(thisdraw);
   draws(thisdraw).rt           = rt;
   
   if abort; fclose('all');break; end 

end % end of draw for loop

% Do we send a trigger at this point? end of sequence trigger:
if EEG == 1 
    io64(ioObj, address, trigger103)
    WaitSecs(triggerdur);
    io64(ioObj, address, 0)
end

% UPDATE BALANCE 
% 1. First subtract the penalty (0.25p for every draw)
balance = balance - (penalty * draw_count);

% 2. Now add winnings or subtract losses
if accuracy == 1
    balance = balance + win; 
elseif accuracy == 0 
    balance = balance - win;    
end

% update the current balance
set.balance         = balance;

% update the position/rate for this trial (if the subject didn't give a
% rate then this should be nan
if answer == 1 || answer == 2
    thisrate        = position;    
else 
    position        = nan;
    thisrate        = position;
    rate_rt         = nan;
end
   
% add trial related info here 
blocktrials.session      = thisession;
blocktrials.block        = thisblock;
blocktrials.trialnumber  = thistrial;
blocktrials.trialonset   = trialstart;
blocktrials.urntype      = thisurn;
blocktrials.sequence     = sequence(1,:);
blocktrials.draws        = draw_count;
blocktrials.response     = answer;
blocktrials.accuracy     = accuracy;
blocktrials.balance      = balance;
blocktrials.thisrate     = thisrate;
blocktrials.ratert       = rate_rt;

WaitSecs(2); %  1= force wait for actual pulse; 0=return this many ms after pulse

% store draws and trials info in the logs mat file 
logs.trialstart         = trialstart;
logs.draws              = draws;
set.blocktrials         = blocktrials;

sub_drawslog             = fullfile(resfolder,sprintf(logs.drawslog,sub,taskname,thisblock,thistrial,thisession));
save(sub_drawslog,'logs');

end