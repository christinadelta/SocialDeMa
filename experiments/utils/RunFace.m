function [set, logs] = RunFace(set, scrn, logs)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
taskNb          = set.taskNb;       % number of task (needed to run the correct task)
blocktrials     = set.blocktrials;  % total number of trials
thisblock       = set.iBlock;       % this is the current block number 
phase           = set.phase;        % 

isi             = set.isi;          % interstimulus interval
jitter          = set.jitter; 
fixation        = set.fixation;     % draw fixation
EEG             = set.EEG;          % is it a EEG session? 

textures        = set.textures;     % this will be used to draw the textures on screen

esckey          = set.code21;       % allows subject to abort experiment;

window          = scrn.window;      % main window
windrect        = scrn.windrect;
xcenter         = scrn.xcenter;
ycenter         = scrn.ycenter;
objectx         = scrn.objectx;
objecty         = scrn.objecty;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;        % slack is ifi/2 (important for timing)
white           = scrn.white;
grey            = scrn.grey;
fixsize         = scrn.fixationsize;
textfont        = scrn.textfont;
textsize        = scrn.textsize;

thisession      = logs.sess;
sub             = logs.sub;
taskname        = logs.task;
resfolder       = logs.resultsfolder;

% create fixation cross offscreen and paste later (faster)
fixationdisplay = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixationdisplay, grey);
Screen('TextFont',fixationdisplay, textfont);
Screen('TextSize',fixationdisplay, fixsize);
DrawFormattedText(fixationdisplay, fixation, 'center', ycenter, white);

% Compute destination rectangle location
destrect        = [xcenter-objectx/4, ycenter-objecty/4, xcenter+objectx/4, ycenter+objecty/4];
set.destrect    = destrect;

%% ----- Run the best-choice economic task ------ %%

abort       = 0;
HideCursor;

if phase == 1
     
    % UNPACK PHASE ONE PARAMETERS 
    sequence        = set.sequence; % this is the current sequence/trial 
    trials          = [];           % store trial info
    
    % START THE BLOCK WITH A FIXATION CROSS
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
    fliptime    = Screen('Flip', window); % flip fixation window
    trialstart  = fliptime;
    
    % object offset
    object_offset   = trialstart + isi + randperm(jitter*1000,1)/1000 - ifi;
    
    % LOOP OVER TRIALS 
    for iTrial = 1: blocktrials
        
        % allow subject to abort the experiment 
        [keyisdown,~,keycode] = KbCheck;
        if keyisdown && keycode(esckey)
            abort = 1;
            break;
        end
        
        % get current trial details 
        thisitem            = sequence(iTrial); % index of the current item 
        
        set.thisitem        = thisitem;
        set.object_offset   = object_offset;
        
        % NOW DARW THE CONTRACT AND SCALE WITHOUT THE SLIDER AND ALLOW THE SUBJECT TO
        % CLICK ONCE TO INIT THE SLIDER
        set = MakeSlider(scrn, set); % first draw the scale and allow subject to click the mouse
        WaitSecs(0.3)                % delay to prevent fast mouse clicks mix 
        
        stimonset           = set.object_onset;

        set = RunSlider(scrn, set);  % once subject click the first time, display slider
        WaitSecs(0.2)

        % UNPACK SETTINGS
        object_offset       = set.object_offset;
        rate_rt             = set.rate_rt;
        position            = set.position;
        
        % BRING FIXATION BACK ON
        Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
        fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     
        object_offset       = fixation_onset + isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
        
        % SAVE TRIAL INFO
        trials(iTrial).sub          = sub;
        trials(iTrial).trialNb      = iTrial;
        trials(iTrial).session      = thisession;
        trials(iTrial).block        = thisblock;
        trials(iTrial).trialstart   = trialstart;
        trials(iTrial).stimonset    = stimonset;
        trials(iTrial).thisitem     = thisitem;
        trials(iTrial).response     = position;
        trials(iTrial).rt           = rate_rt;
        
        if abort; fclose('all');break; end 
        
    end % end of trials loop
    
    WaitSecs(1); % wait two sec before flipping to the next block/

    logs.trials              = trials; % save the samples 

    sublogs                  = fullfile(resfolder,sprintf(logs.trialog,sub,taskname,thisblock,thisession,phase));
    save(sublogs,'logs');
    
else % if this is phase 2 
    
    % UNPACK SCREEN STRUCT for phase 2
    textbold        = scrn.textbold;
    smalltext       = scrn.smalltext;
    baserect        = scrn.baserect;
    centeredrect    = scrn.centeredrect;
    penwidth        = scrn.penwidth;
    orange          = scrn.orange;
    
    % CREATE RECTANGLES
    % BASE RECTANGLES -- white rect
    whiterect_window = Screen('OpenOffscreenWindow',window);
    Screen('TextSize', whiterect_window, textsize);
    Screen('FillRect', whiterect_window, grey ,windrect);
    Screen('FrameRect', whiterect_window, white, centeredrect, penwidth)

    % BASE RECTANGLES -- orange rect
    orangerect_window = Screen('OpenOffscreenWindow',window);
    Screen('FillRect', orangerect_window, grey ,windrect);
    Screen('FrameRect', orangerect_window, orange, centeredrect, penwidth);
    
    % UNPACK TRIGGER STUFF
    if EEG == 1 
    
        ioObj       = set.ioObject;
        status      = set.status;
        triggerdur  = set.triggerdur;
        address     = set.address;
        
%         trigger11   = set.trigger11;
%         trigger12   = set.trigger12;

        trigger1    = set.trigger1;
        trigger13   = set.trigger13;
        trigger14   = set.trigger14;
        trigger15   = set.triiger15;
        
        trigger100  = set.trigger100; % start of sequence
        trigger101  = set.trigger101; % end of sequence
        
    end
    
    % UNPACK SETTINGS STRUCT for phase 2
    thistrial       = set.thisTrial;    % current sequence/trial number
    sequence        = set.sequence;     % this is the current sequence/trial
    fixduration     = set.fix_dur;      % fixation duration 
    response        = set.response;     % response time (5 sec or self-paced)
    feedbacktime    = set.feedback;     % feedback duration
    smalltex        = set.smalltex;     % small textures 
    stimduration    = set.stimdur;      % stimulus duration

     % UNPACK RESPONSE KEYS
    code1           = set.code1;
    code2           = set.code2;

    samples         = set.samples;
    trials          = [];   % store sequence info
    previous        = [];   % store previous samples
    blocktrials     = [];   % here we store the info of the current trial
    xcntr           = xcenter - 60; % works better than 'center'
        
    % START THE TRIAL WITH A FIXATION CROSS
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
    fliptime    = Screen('Flip', window); % flip fixation window
    trialstart  = fliptime;
    
    % send sequence start trigger
    if EEG == 1 
        io64(ioObj, address, trigger100)
        WaitSecs(triggerdur);
        io64(ioObj, address, 0) % return port to zero
    end
    
    % object offset
    object_offset   = trialstart + isi + randperm(jitter*1000,1)/1000 - ifi;
    
    % start sampling 
    for s = 1:samples
        
        % first allow subject to abort the experiment (if they press esc)
        [keyisdown,~,keycode] = KbCheck;
        if keyisdown && keycode(esckey)
            abort = 1;
            break;
        end
        
        thisitem        = sequence(1,s);    % index of the current item       
        Xchange         = xcntr-60;         % this will help define x positions of the small images
        Ychange         = ycenter;          % this is the y position of the small images
        xs              = Xchange;
        ys              = Ychange;
        tmp             = 0;
        
        % DISPLAY CONTRACT 
        % FIRST DRAW THE RECT
        Screen('CopyWindow', whiterect_window,window, windrect, windrect)
        Screen('TextSize', window, textsize);
        Screen('TextStyle', window, 0)
        
        if s == 1
            % DISPLAY THE CURRENT FACE - Backround window (the number of
            % sample)
            DrawFormattedText(window,sprintf('Face %d/10',s), 'center', ycenter-200, white);

            % DISPLAY THE CURRENT FACE - the actual face
            Screen('DrawTexture',window, textures{thisitem}, [], destrect);     % display thisitem
            object_onset    = Screen('Flip', window, object_offset - slack);            % here the current image (thisitem) is fliped
            
        else % if it is not the first sample
            
            previous_len = length(previous);  % how many previous faces?
            
            % image size and position on screen (for the small images)
            smallrects  = nan(4,previous_len); 
            
            for l = 1:previous_len
                
                indx            = previous_len - tmp; % this will be used to present the last bead (of the previous beads) first
                smallrects(:,indx) = [Xchange-objectx/9, Ychange-objecty/9, Xchange+objectx/9, Ychange+objecty/9];
                
                Xchange = Xchange - 60;
                tmp         = tmp + 1;  % update tmp 
            end
            
            % DISPLAY THE CURRENT FACE - Backround window (the number of
            % sample)
            DrawFormattedText(window,sprintf('Face %d/10',s), 'center', ycenter-200, white);
           
            % DISPLAY THE CURRENT FACE - the actual face
            Screen('DrawTexture', window, textures{thisitem}, [], destrect);     % display thisitem
            Screen('DrawTextures', window, [smalltex{previous}], [], smallrects);     % display small items
            object_onset    = Screen('Flip', window, object_offset - slack);            % here the current image (thisitem) is fliped
        end
    
        % send stimulus trigger
        if EEG == 1 
            io64(ioObj, address, trigger1)
            WaitSecs(triggerdur);
            io64(ioObj, address, 0) % return port to zero
        end
        
        object_offset   = object_onset + stimduration - ifi; % add jitter here?
        
        tmp             = 0; % update tmp for the coloured rect now
        % 2. SHOW GO SIGNAL
        % SHOW RECT WITH THE LIST OF BEADS(S)
        Screen('CopyWindow', orangerect_window,window, windrect, windrect)
        Screen('TextSize', window, textsize);
        Screen('TextStyle', window, 0) % NOT BOLD
        
        DrawFormattedText(window,sprintf('Face %d/10',s), 'center', ycenter-200, white);
        DrawFormattedText(window, '1:', xcntr-180, ycenter+220, orange);
        DrawFormattedText(window, 'accept date?', xcntr-230, ycenter+250, orange);
        DrawFormattedText(window, '2:', xcntr+200, ycenter+220, orange);
        DrawFormattedText(window, 'sample again?', xcntr+140, ycenter+250, orange);
        % DISPLAY THE CURRENT FACE - the actual face
        Screen('DrawTexture',window, textures{thisitem}, [], destrect);     % display thisitem
        
        if s > 1
            
            % image size and position on screen (for the small images)
            smallrects  = nan(4,previous_len);
            
            for l = 1:previous_len
                indx            = previous_len - tmp; % this will be used to present the last bead (of the previous beads) first
                smallrects(:,indx) = [xs-objectx/9, ys-objecty/9, xs+objectx/9, ys+objecty/9];

                xs      = xs-60;
                tmp     = tmp + 1;
            end
            
            % draw the previous faces
            Screen('DrawTextures', window, [smalltex{previous}], [], smallrects);     % display small items   
        end
        
        % flip screen
        response_onset    = Screen('Flip', window, object_offset - slack);
        
        % send response prompt trigger
        if EEG == 1 
            io64(ioObj, address, trigger15)
            WaitSecs(triggerdur);
            io64(ioObj, address, 0) % return port to zero
        end
        
        % WAIT FOR RESPONSE 
        rt                  = NaN;
        answer              = NaN;
        resp_input          = 0;
        responseTrigNotSent = 1;
        
        while resp_input == 0 && (GetSecs - response_onset) < response - 2*slack
            
            [~, secs, keycode] = KbCheck; % check for input
            
            if keycode(1,code1) % 
                resp_input  = code1;
                rt          = secs - response_onset;
                answer      = 1; % subject accepted a face/date
                respmade    = secs;
                
                % send response trigger -- subject accepted an option
%                 if EEG == 1 && responseTrigNotSent == 1
%                     io64(ioObj, address, trigger11)
%                     WaitSecs(triggerdur);
%                     io64(ioObj, address, 0) % return port to zero
%                     responseTrigNotSent = 0;
%                 end
                
            elseif keycode(1,code2) %  
                resp_input  = code2;
                rt          = secs - response_onset;
                answer      = 2; % subject sampled again
                respmade    = secs;
                
                % send response trigger -- subject sampled again
%                 if EEG == 1 && responseTrigNotSent==1
%                     io64(ioObj, address, trigger12)
%                     WaitSecs(triggerdur);
%                     io64(ioObj, address, 0) % return port to zero 
%                     responseTrigNotSent = 0;
%                 end
%                 
            else
                resp_input  = 0; 
                rt          = nan; 
                answer      = nan; % no response 
                respmade    = GetSecs;
            end
        end % end of response while loop
        
        % object offset 
        object_offset   = respmade + isi - ifi;                             % contract window self paced or on for 5000 ms
       
        % IF SUBJECT ACCEPTED A FACE/DATE SHOW FFEDBACK AND BREAK SEQUENCE
        % IF SUBJECT CHOSE TO SAMPLE AGAIN, SHOW FIXATION AND MOVE TO THE
        % NEXT FACE
        if answer == 1
            
            % BRING FIXATION BACK ON
            Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
            fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

            object_offset       = fixation_onset + fixduration - ifi;  % add jitter here?
            
            % DISPLAY CHOSEN CONTRACT
            feedback_window = Screen('OpenOffscreenWindow', window, windrect);
            Screen('TextSize', feedback_window, textsize);
            Screen('FillRect', feedback_window, grey ,windrect);
            DrawFormattedText(feedback_window, 'Congratulations! This is your new date.', 'center', ycenter-250, white);

            % DISPLAY THE CURRENT FACE
            Screen('DrawTexture', feedback_window, textures{thisitem}, [], destrect);     % display the chosen item
            Screen('CopyWindow',feedback_window, window, windrect, windrect);
            object_onset        = Screen('Flip', window, object_offset - slack);    % flip window
            
            % send feedback screen trigger (accepted date trigger)
            if EEG == 1
                io64(ioObj, address, trigger13)
                WaitSecs(triggerdur);
                io64(ioObj, address, 0) % return port to zero
            end
            
            % DISPLAY OFFSET 
            object_offset   = object_onset + feedbacktime - ifi;
            
            % BRING FIXATION BACK ON  
            Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
            fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

            object_offset   = fixation_onset + fixduration + isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
            
            break; % break from sequence 
            
        else % if subject wants to sample again or if subject doesn't give a response, move to the next sample
            
            % if by the 10th sample subject chooses to sample again 
            if s == samples 
                
                % BRING FIXATION BACK ON
                Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
                fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

                object_offset       = fixation_onset + fixduration - ifi; % add jitter here?

                feedback_sampling = Screen('OpenOffscreenWindow',window);
                Screen('TextSize', feedback_sampling, textsize);
                Screen('FillRect', feedback_sampling, grey ,windrect);
                DrawFormattedText(feedback_sampling, 'Oh No! :(', 'center', ycenter-50, white);
                DrawFormattedText(feedback_sampling, 'You you are not allowed to sample again', 'center', ycenter, white);
                Screen('CopyWindow',feedback_sampling, window, windrect, windrect);
                object_onset        = Screen('Flip', window, object_offset - slack);    % flip window
                
                % send feedback screen trigger - you can't sample again!
                if EEG == 1
                    io64(ioObj, address, trigger14)
                    WaitSecs(triggerdur);
                    io64(ioObj, address, 0) % return port to zero
                end
            
                % DISPLAY REQUESTED OBJECT OFFSET 
                object_offset       = object_onset + feedbacktime - ifi;
                
                % BRING FIXATION BACK ON
                Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
                fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

                object_offset       = fixation_onset + fixduration - ifi; % add jitter here?
                
                % DISPLAY CHOSEN FACE/DATE
                feedback_window = Screen('OpenOffscreenWindow', window, windrect);
                Screen('TextSize', feedback_window, textsize);
                Screen('FillRect', feedback_window, grey ,windrect);
                DrawFormattedText(feedback_window, 'Congratulations! This is your new date.', 'center', ycenter-250, white);

                % DISPLAY THE CURRENT FACE
                Screen('DrawTexture', feedback_window, textures{thisitem}, [], destrect);     % display the chosen item
                Screen('CopyWindow',feedback_window, window, windrect, windrect);
                object_onset        = Screen('Flip', window, object_offset - slack);    % flip window

                % send feedback screen trigger
                if EEG == 1
                    io64(ioObj, address, trigger13)
                    WaitSecs(triggerdur);
                    io64(ioObj, address, 0) % return port to zero
                end
                
                % DISPLAY OBJECT OFFSET 
                object_offset   = object_onset + feedbacktime - ifi;
                
            end % end of s statement  
            
            % BRING FIXATION BACK ON  
            Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
            fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

            fprintf('prompt was on for %3.4f\n', fixation_onset - object_onset);        % time interval from the the flip of the contract until fixation

            object_offset   = fixation_onset + fixduration + isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
            
        end % end of feedback statement
        
        % store the cuurent face to show at the bottof of the screen during
        % next samples
        previous(s)            = thisitem;
        
        % save the sequence-sampling info 
        trials(s).session      = thisession;
        trials(s).block        = thisblock;
        trials(s).trialnumber  = thistrial;
        trials(s).trialonset   = trialstart;
        trials(s).thisitem     = thisitem;
        trials(s).rt           = rt;
        
        if abort; fclose('all');break; end 
    end % end of sampling loop
    
    % store the sequence info in logs 
    logs.trials                 = trials; % save the samples 

    sublogs                     = fullfile(resfolder,sprintf(logs.trialog,sub,taskname,thisblock,thistrial,thisession,phase));
    save(sublogs,'logs');
    
    % send sequence end trigger
    if EEG == 1 
        io64(ioObj, address, trigger101)
        WaitSecs(triggerdur);
        io64(ioObj, address, 0) % return port to zero
    end

    % save the current trial info
    blocktrials.session         = thisession;
    blocktrials.block           = thisblock;
    blocktrials.trialnumber     = thistrial;
    blocktrials.trialonset      = trialstart;
    blocktrials.sequence        = sequence;
    blocktrials.numsamples      = s;
    blocktrials.chosenitem      = thisitem; 

    set.blocktrials             = blocktrials; % (save the info of this trial) this will go to the main script
    
    WaitSecs(1); % wait two sec before flipping to the next block/

end % end of phase statement 


end