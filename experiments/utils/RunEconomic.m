function [set, logs] = RunEconomic(set, scrn, logs)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
blocktrials     = set.blocktrials;  % total number of trials
thisblock       = set.iBlock;    % this is the current block number 
phase           = set.phase;        % 

isi             = set.isi;          % interstimulus interval
jitter          = set.jitter; 
fixation        = set.fixation;     % draw fixation
EEG             = set.EEG;          % is it a EEG session? 

esckey          = set.code21;       % allows subject to abort experiment;

window          = scrn.window;      % main window
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

%% ----- Run the best-choice economic task ------ %%

abort       = 0;
HideCursor;

if phase == 1
     %% RUN PHASE TWO %%
    
    % UNPACK PHASE ONE PARAMETERS 
    sequence        = set.sequence;     % this is the current sequence/trial
    fixduration     = set.fix_dur;       % fixation duration 
   
    % UNPACK STIMULI 
    price           = set.price;

    trials          = [];   % store trial info
 
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
        thisprice           = price(thisitem);
        
        % Update settings struct
        set.thisprice       = thisprice;
        set.object_offset   = object_offset;
 
        % NOW DARW THE CONTRACT AND SCALE WITHOUT THE SLIDER AND ALLOW THE SUBJECT TO
        % CLICK ONCE TO INIT THE SLIDER
        set                 = MakeSlider(scrn, set); % first draw the scale and allow subject to click the mouse
        WaitSecs(0.3)                % delay to prevent fast mouse clicks mix 
        
        stimonset           = set.object_onset;

        set                 = RunSlider(scrn, set);  % once subject click the first time, display slider
        WaitSecs(0.2)

        % UNPACK SETTINGS
        object_offset       = set.object_offset;
        rate_rt             = set.rate_rt;
        position            = set.position;
       
        % PUT FIXATION BACK ON
        Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
        fixation_onset  = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     
        object_offset   = fixation_onset + isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
        
        % SAVE TRIAL INFO
        trials(iTrial).sub          = sub;
        trials(iTrial).trialNb      = iTrial;
        trials(iTrial).session      = thisession;
        trials(iTrial).block        = thisblock;
        trials(iTrial).trialstart   = stimonset;
        trials(iTrial).thisitem     = thisitem;
        trials(iTrial).thisprice    = thisprice;
        trials(iTrial).response     = position;
        trials(iTrial).rt           = rate_rt;
        
        if abort; fclose('all');break; end 
        
    end % end of trial loop
    
    logs.trials              = trials; % save the sequence

    sublogs                  = fullfile(resfolder,sprintf(logs.trialog,sub,taskname,thisblock,thisession,phase));
    save(sublogs,'logs');
    
    WaitSecs(1); % wait two sec before flipping to the next block
    

else % if phase == 2
    %% RUN PHASE TWO %%
      
    % UNPACK SETTINGS STRUCT for phase 2
    thistrial       = set.thisTrial;    % current sequence/trial number
    sequence        = set.sequence;     % this is the current sequence/trial
    fixduration     = set.fix_dur;       % fixation duration 
    response        = set.response;     % response time (5 sec or self-paced)
    feedbacktime    = set.feedback;     % feedback duration
    stimduration    = set.stimdur;      % stimulus duration
    balance         = set.balance;      % current balance
    rewards         = set.reward;
    
    textbold        = scrn.textbold;
    smalltext       = scrn.smalltext;
    baserect        = scrn.baserect;
    centeredrect    = scrn.centeredrect;
    penwidth        = scrn.penwidth;
    orange          = scrn.orange;
    
    % UNPACK STIMULI 
    prices          = set.price;

     % UNPACK RESPONSE KEYS
    code1           = set.code1; % accept price 
    code2           = set.code2; % sample again

    samples         = set.samples; % total number of samples
        
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
        
        trigger11   = set.trigger11;
        trigger12   = set.trigger12;
        trigger13   = set.trigger13;
        trigger14   = set.trigger14;
        trigger15   = set.triiger15;
        
        trigger100  = set.trigger100; % start of sequence
        trigger101  = set.trigger101; % end of sequence
        
    end
    
    trial_samples   = [];     % store trial info
    previous        = [];     % store the previous prices here to show at the bottom of the screen
    blocktrials     = [];     % here we store the info of the current sequence
    sequenceprices  = nan(1,samples);
    xcntr           = xcenter - 80; % works better than 'center'
    
    % convert sequence prices to num and save them to get the ranks
    for p = 1:samples
        temp                    = sequence(p);
        sequenceprices(p)       = prices(temp);
    end
    
    % GET RANKS (three lowest prices)
    minprices                   = mink(sequenceprices,3);
 
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
        
        if EEG == 1
            trigger1 = 0 + s; % contract/price trigger codes [1:10]
        end
        
        thisitem        = sequence(s); % index of the current item
        thisprice       = prices(thisitem);
        xcenters        = xcntr - 60; % start adding the previous prices at xcenter - 60
        xs              = xcenters;
        tmp             = 0;
        
        % DISPLAY CONTRACT 
        % FIRST DRAW THE RECT
        Screen('CopyWindow', whiterect_window,window, windrect, windrect)
        Screen('TextSize', window, textsize);
        
        if s == 1
            
            % if this is the 1st sample show only the 1st contract
            DrawFormattedText(window,sprintf('Contract %d/10',s), 'center', ycenter-200, white);
            Screen('TextStyle', window, textbold);
            DrawFormattedText(window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter, white); 

            object_onset        = Screen('Flip', window, object_offset - slack);    % flip window

        else
            
            previous_len        = length(previous); % how many previous prices?
            
            % background_window = Screen('OpenOffscreenWindow', window, windrect);
            DrawFormattedText(window, sprintf('Contract %d/10',s), 'center', ycenter-200, white);
            Screen('TextStyle', window, textbold);
            DrawFormattedText(window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter, white); 
            
            for l = 1:previous_len
                
                indx = previous_len - tmp; % this will be used to present the last bead (of the previous beads) first
                % show the previous prices at the bottom of the screen
                Screen('TextSize', window, smalltext);
                Screen('TextStyle', window, 0)
                DrawFormattedText(window, sprintf('£%3.2f', previous{indx}), xcenters, ycenter, white); 
                
                % update xcenters, so that previous prices are not
                % displayed on top of each other
                xcenters = xcenters - 60;
                tmp = tmp + 1;
                
            end  
            object_onset = Screen('Flip', window, object_offset - slack);    % flip window
            
        end
        
        % send stimulus trigger
        if EEG == 1 
            io64(ioObj, address, trigger1)
            WaitSecs(triggerdur);
            io64(ioObj, address, 0) % return port to zero
        end
        
        object_offset   = object_onset + stimduration - ifi; % add jitter here?
        
        tmp             = 0; % update tmp for the colour now
        
        % 2. SHOW GO SIGNAL
        % SHOW RECT WITH THE LIST OF BEADS(S)
        Screen('CopyWindow', orangerect_window,window, windrect, windrect)
        Screen('TextSize', window, textsize);
        Screen('TextStyle', window, 0) % NOT BOLD
        DrawFormattedText(window,sprintf('Contract %d/10',s), 'center', ycenter-200, white);
        DrawFormattedText(window, '1:', xcntr-180, ycenter+220, orange);
        DrawFormattedText(window, 'accept option?', xcntr-230, ycenter+250, orange);
        DrawFormattedText(window, '2:', xcntr+200, ycenter+220, orange);
        DrawFormattedText(window, 'sample again?', xcntr+140, ycenter+250, orange);
        Screen('TextSize', window, smalltext);
        Screen('TextStyle', window, 0) % NOT BOLD
        DrawFormattedText(window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter, white); 
        
        if s > 1
            
            previous_len = length(previous);
            
            for l = 1:previous_len
                
                indx = previous_len - tmp; % this will be used to present the last bead (of the previous beads) first
                % show the previous prices 
                DrawFormattedText(window, sprintf('£%3.2f', previous{indx}), xs, ycenter, white); 
                % update xcenter, so that previous prices are not
                % displayed on top of each other
                xs = xs - 60;
                tmp     = tmp + 1;
            end
        end

        response_onset    = Screen('Flip', window, object_offset - slack);
        
        % send sequence start trigger
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
                answer      = 1; % subject accepted a contract
                respmade    = secs;
                
                % send response trigger -- subject accepted an option
                if EEG == 1 && responseTrigNotSent==1
                    io64(ioObj, address, trigger11)
                    WaitSecs(triggerdur);
                    io64(ioObj, address, 0) % return port to zero   
                end
                
            elseif keycode(1,code2) %  
                resp_input  = code2;
                rt          = secs - response_onset;
                answer      = 2; % subject sampled again
                respmade    = secs;
                
                % send response trigger -- subject sampled again
                if EEG == 1 && responseTrigNotSent==1
                    io64(ioObj, address, trigger12)
                    WaitSecs(triggerdur);
                    io64(ioObj, address, 0) % return port to zero   
                end
               
            else
                resp_input  = 0; 
                rt          = nan; 
                answer      = nan; % no response 
                respmade    = GetSecs;
            end
        end % end of response while loop
        
        % object offset 
        object_offset   = respmade + isi - ifi;                             % contract window self paced or on for 5000 ms
       
         % what is the reward and rank of the accepted?
        numprice        = thisprice;
        
        if numprice == minprices(1)
            thisreward  = rewards(1);
            thisrank    = 1;
            
        elseif numprice == minprices(2)
            thisreward  = rewards(2);
            thisrank    = 2;
            
        elseif numprice == minprices(3)
            thisreward  = rewards(3);
            thisrank    = 3;
            
        else
            thisreward  = 0; 
            thisrank    = 0;
        end
        
        % IF SUBJECT ACCEPTED A CONTRACT SHOW FFEDBACK AND BREAK SEQUENCE
        % IF SUBJECT CHOSE TO SAMPLE AGAIN, SHOW FIXATION AND MOVE TO THE
        % NEXT ITEM
        if answer == 1
            
            % BRING FIXATION BACK ON
            Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
            fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

            object_offset   = fixation_onset + fixduration - ifi; % add jitter here?
            
            % DISPLAY FEEDBACK: SHOW THE ACCEPTED CONTRACT AND REWARD 
            if thisrank == 0
 
                % DISPLAY CHOSEN CONTRACT
                feedback_window = Screen('OpenOffscreenWindow', window, windrect);
                Screen('TextSize', feedback_window, textsize)
                Screen('TextStyle', window, 0) % NOT BOLD;
                Screen('FillRect', feedback_window, grey ,windrect);
                DrawFormattedText(feedback_window, 'Congratulations! This is the price  of your new smartphone contract.', 'center', ycenter-200, white);
                DrawFormattedText(feedback_window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter, white); 
            
            else
                % DISPLAY CHOSEN CONTRACT
                feedback_window = Screen('OpenOffscreenWindow', window, windrect);
                Screen('TextSize', feedback_window, textsize);
                Screen('TextStyle', window, 0) % NOT BOLD
                Screen('FillRect', feedback_window, grey ,windrect);
                DrawFormattedText(feedback_window, 'Congratulations! This is the price  of your new smartphone contract.', 'center', ycenter-200, white);
                DrawFormattedText(feedback_window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter-50, white); 
                DrawFormattedText(feedback_window, sprintf('Your reward is %3.3f credits',thisreward), 'center', scrn.ycenter, scrn.white);
            end
            
            Screen('CopyWindow',feedback_window, window, windrect, windrect);
            object_onset        = Screen('Flip', window, object_offset - slack);    % flip window
            
            % send feedback screen trigger
            if EEG == 1
                io64(ioObj, address, trigger13)
                WaitSecs(triggerdur);
                io64(ioObj, address, 0) % return port to zero
            end
            
            % DISPLAY REQUESTED OBJECT OFFSET 
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

                object_offset   = fixation_onset + fixduration - ifi; % add jitter here?

                % DISPLAY FEEDBACK (WARNING) 
                feedback_sampling = Screen('OpenOffscreenWindow',window);
                Screen('TextSize', feedback_sampling, textsize);
                Screen('TextStyle', window, 0) % NOT BOLD
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
                object_offset   = object_onset + feedbacktime - ifi;
                
                % BRING FIXATION BACK ON
                Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
                fixation_onset  = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

                object_offset   = fixation_onset + fixduration - ifi; % add jitter here?
                
                % DISPLAY FEEDBACK: SHOW THE ACCEPTED CONTRACT AND REWARD 
                if thisrank == 0

                    % DISPLAY CHOSEN CONTRACT
                    feedback_window = Screen('OpenOffscreenWindow', window, windrect);
                    Screen('TextSize', feedback_window, textsize);
                    Screen('TextStyle', window, 0); % NOT BOLD
                    Screen('FillRect', feedback_window, grey ,windrect);
                    DrawFormattedText(feedback_window, 'Congratulations! This is the price  of your new smartphone contract.', 'center', ycenter-200, white);
                    DrawFormattedText(feedback_window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter, white); 

                else
                    % DISPLAY CHOSEN CONTRACT
                    feedback_window = Screen('OpenOffscreenWindow', window, windrect);
                    Screen('TextSize', feedback_window, textsize);
                    Screen('TextStyle', window, 0) % NOT BOLD
                    Screen('FillRect', feedback_window, grey ,windrect);
                    DrawFormattedText(feedback_window, 'Congratulations! This is the price  of your new smartphone contract.', 'center', ycenter-200, white);
                    DrawFormattedText(feedback_window, sprintf('Price: £%3.2f', thisprice), 'center', ycenter-50, white); 
                    DrawFormattedText(feedback_window, sprintf('Your reward is %3.3f credits',thisreward), 'center', scrn.ycenter, scrn.white);
                end

                Screen('CopyWindow',feedback_window, window, windrect, windrect);
                object_onset     = Screen('Flip', window, object_offset - slack); % flip window
                
                % send feedback screen trigger
                if EEG == 1
                    io64(ioObj, address, trigger13)
                    WaitSecs(triggerdur);
                    io64(ioObj, address, 0) % return port to zero
                end

                % DISPLAY OFFSET 
                object_offset    = object_onset + feedbacktime - ifi;
                
            end
            
            % BRING FIXATION BACK ON AND MOVE TO THE NEXT TRIAL/SAMPLE
            Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
            fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     
            object_offset       = fixation_onset + fixduration + isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
            
        end
        
        % store the current price to show at the bottom of the screen
        % (during the next samples)
        previous{s}                   = thisprice;
        
        % save the sequence-sampling info 
        trial_samples(s).session      = thisession;
        trial_samples(s).block        = thisblock;
        trial_samples(s).trialnumber  = thistrial;
        trial_samples(s).trialonset   = trialstart;
        trial_samples(s).thisitem     = thisitem;
        trial_samples(s).thisprice    = thisprice;
        trial_samples(s).rt           = rt;
        
        if abort; fclose('all');break; end 
        
    end % end of sampling for loop  
    
    logs.trialsamples        = trial_samples; % save the samples 

    sublogs                  = fullfile(resfolder,sprintf(logs.trialog,sub,taskname,thisblock,thistrial,thisession,phase));
    save(sublogs,'logs');
    
    % send sequence end trigger
    if EEG == 1 
        io64(ioObj, address, trigger101)
        WaitSecs(triggerdur);
        io64(ioObj, address, 0) % return port to zero
    end
    
    % update balance 
    balance                  = balance + thisreward;
    set.balance              = balance;
    chosenitem               = thisitem;

    % save the current trial info
    blocktrials.session      = thisession;
    blocktrials.block        = thisblock;
    blocktrials.trialnumber  = thistrial;
    blocktrials.trialonset   = trialstart;
    blocktrials.sequence     = sequence';
    blocktrials.numsamples   = s;
    blocktrials.chosenitem   = chosenitem;
    blocktrials.chosenprice  = thisprice;
    blocktrials.rank         = thisrank;
    blocktrials.reward       = thisreward;
    blocktrials.balance      = balance;

    set.blocktrials          = blocktrials; % (save the info of this trial) this will go to the main script
    
    WaitSecs(1); % wait two sec before flipping to the next block/

    
end % end of phase statement


end