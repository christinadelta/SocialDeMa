function [set, logs] = RunFace(set, scrn, logs)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
taskNb          = set.taskNb;       % number of task (needed to run the correct task)
blocktrials     = set.blocktrials;  % total number of trials
thisblock       = set.iBlock;    % this is the current block number 
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
destrect        = [xcenter-objectx/2, ycenter-objecty/2, xcenter+objectx/2, ycenter+objecty/2];
%destrect        = [xcenter-objectx/4, ycenter-objecty/4, xcenter+objectx/4, ycenter+objecty/4];


%% ----- Run the best-choice economic task ------ %%

abort       = 0;
HideCursor;

if phase == 1
     
    % UNPACK PHASE ONE PARAMETERS 
    sequence        = set.sequence;     % this is the current sequence/trial
    fixduration     = set.fixdur;       % fixation duration 
    response        = set.response;     % response time (5 sec or self-paced)
    
     % UNPACK RESPONSE KEYS
    code1       = set.code1;
    code2       = set.code2;
    code3       = set.code3;
    code4       = set.code4;
    code5       = set.code5;
    code6       = set.code6;
    code7       = set.code7;
    code8       = set.code8;
    code9       = set.code9;
    
    trials      = [];   % store trial info
    
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
        thisitem        = sequence(iTrial); % index of the current item 
        
        % DISPLAY THE CURRENT FACE
        Screen('DrawTexture', window, textures{thisitem}, [], destrect);     % display thisitem
        object_onset    = Screen('Flip', window, object_offset - slack);            % here the current image (thisitem) is fliped
        
        rt                  = NaN;
        answer              = NaN;
        resp_input          = 0;
        
        while resp_input == 0 && (GetSecs - object_onset) < response - 2*slack
            
            [~, secs, keycode] = KbCheck; % check for input
            
            if keycode(1,code1) % if subject chose the green urn 
                resp_input  = code1;
                rt          = secs - object_onset;
                answer      = 1; % green urn
                respmade    = secs;
                
            elseif keycode(1,code2) % if subject chose the green urn 
                resp_input  = code2;
                rt          = secs - object_onset;
                answer      = 2; % green urn
                respmade    = secs;
                
            elseif keycode(1,code3) % if subject chose to draw again
                resp_input  = code3;
                rt          = secs - object_onset;
                answer      = 3; % draw-again 
                respmade    = secs;
                
            elseif keycode(1,code4)
                resp_input  = code4;
                rt          = secs - object_onset;
                answer      = 4; % draw-again 
                respmade    = secs;
                
            elseif keycode(1,code5)
                resp_input  = code5;
                rt          = secs - object_onset;
                answer      = 5; % draw-again 
                respmade    = secs;
                
            elseif keycode(1,code6)
                resp_input  = code6;
                rt          = secs - object_onset;
                answer      = 6; % draw-again 
                respmade    = secs;
                
            elseif keycode(1,code7)
                resp_input  = code7;
                rt          = secs - object_onset;
                answer      = 7; % draw-again 
                respmade    = secs;
                
            elseif keycode(1,code8)
                resp_input  = code8;
                rt          = secs - object_onset;
                answer      = 8; % draw-again 
                respmade    = secs;
                
            elseif keycode(1,code9)
                resp_input  = code9;
                rt          = secs - object_onset;
                answer      = 9; % draw-again 
                respmade    = secs;
                
            else
                resp_input = 0; 
                rt          = nan; 
                answer      = nan; % answered A
                respmade    = GetSecs;
            end
        end % end of response while loop
        
        % object offset 
        object_offset       = respmade + isi - ifi;                             % contract window self paced or on for 5000 ms
        
        Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
        fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     

        fprintf('prompt was on for %3.4f\n', fixation_onset - object_onset);        % time interval from the the flip of the contract until fixation

        object_offset       = fixation_onset + fixduration + isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
        
        % SAVE TRIAL INFO
        trials(iTrial).sub          = sub;
        trials(iTrial).trialNb      = iTrial;
        trials(iTrial).session      = thisession;
        trials(iTrial).block        = thisblock;
        trials(iTrial).trialstart   = trialstart;
        trials(iTrial).trialstart   = trialstart;
        trials(iTrial).thisitem     = thisitem;
        trials(iTrial).response     = answer;
        trials(iTrial).rt           = rt;
        
        if abort; fclose('all');break; end 
        
    end % end of trials loop
    
    WaitSecs(1); % wait two sec before flipping to the next block/

    logs.trials              = trials; % save the samples 

    sublogs                  = fullfile(resfolder,sprintf(logs.trialog,sub,taskname,thisblock,thisession,phase));
    save(sublogs,'logs');
    
else % if this is phase 2 

     if EEG == 1 
        
        % UNPACK TRIGGERS
        sp          = set.sp;
        trigger11   = set.trigger11;
        trigger12   = set.trigger12;
        trigger13   = set.trigger13;
        
        trigger100  = set.trigger100; % start of sequence
        trigger101  = set.trigger101; % end of sequence

     end
    
      % UNPACK SETTINGS STRUCT for phase 2
    thistrial       = set.thisTrial;    % current sequence/trial number
    sequence        = set.sequence;     % this is the current sequence/trial
    fixduration     = set.fixdur;       % fixation duration 
    response        = set.response;     % response time (5 sec or self-paced)
    feedbacktime    = set.feedback;     % feedback duration
    smalltex        = set.smalltex;     % small textures 
    
    
     % UNPACK RESPONSE KEYS
    code1           = set.code1;
    code2           = set.code2;

    samples         = set.samples;
    trials          = [];   % store trial info
    previous        = nan;   % store previous samples
    blocktrials     = [];   % here we store the info of the current sequence
    
    % START THE TRIAL WITH A FIXATION CROSS
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
    fliptime    = Screen('Flip', window); % flip fixation window
    trialstart  = fliptime;
    
    % send sequence start trigger
    if EEG == 1 
        sp.sendTrigger(trigger100)
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
            trigger1 = 0 + s; % contract triggers 1:10
        end
        
        thisitem        = sequence(s); % index of the current item
        
        Xchange         = xcenter-600; % this will help define x positions of the small images
        Ychange         = ycenter+350; % this is the y position of the small images
        
        if s == 1
            % DISPLAY THE CURRENT FACE - Backround window (the number of
            % sample)
            background_window = Screen('OpenOffscreenWindow', window, windrect);
            Screen('TextSize', background_window, textsize);
            Screen('FillRect', background_window, grey ,windrect);
            DrawFormattedText(background_window,sprintf('Contract %d/10',s), 'center', ycenter-300, white);

             % DISPLAY THE CURRENT FACE - the actual face
            Screen('DrawTexture', background_window, textures{thisitem}, [], destrect);     % display thisitem
            Screen('CopyWindow',background_window, window, windrect, windrect);
            object_onset    = Screen('Flip', window, object_offset - slack);            % here the current image (thisitem) is fliped
            
        else % if it is not the first sample
            
            previous_len = length(previous);  % how many previous faces?
            
            % image size and position on screen (for the small images)
            smallrects  = nan(4,previous_len); 
            
            for l = 1:previous_len
                smallrects(:,l) = [Xchange-objectx/6, Ychange-objecty/6, Xchange+objectx/6, Ychange+objecty/6];
                
                Xchange = Xchange+100;
            end
            
            % DISPLAY THE CURRENT FACE - Backround window (the number of
            % sample)
            background_window = Screen('OpenOffscreenWindow', window, windrect);
            Screen('TextSize', background_window, textsize);
            Screen('FillRect', background_window, grey ,windrect);
            DrawFormattedText(background_window,sprintf('Contract %d/10',s), 'center', ycenter-300, white);

             % DISPLAY THE CURRENT FACE - the actual face
            Screen('DrawTexture', background_window, textures{thisitem}, [], destrect);     % display thisitem
            Screen('DrawTextures', background_window, [smalltex{previous}], [], smallrects);     % display small items
            
            Screen('CopyWindow',background_window, window, windrect, windrect);
            object_onset    = Screen('Flip', window, object_offset - slack);            % here the current image (thisitem) is fliped
        end
    
        % send sequence start trigger
        if EEG == 1
            sp.sendTrigger(trigger1) % blue urn -- high prob blue bead trigger
        end
        
        % WAIT FOR RESPONSE 
        rt                  = NaN;
        answer              = NaN;
        resp_input          = 0;
        responseTrigNotSent = 1;
        
        while resp_input == 0 && (GetSecs - object_onset) < response - 2*slack
            
            [~, secs, keycode] = KbCheck; % check for input
            
            if keycode(1,code1) % 
                resp_input  = code1;
                rt          = secs - object_onset;
                answer      = 1; % subject accepted a face/date
                respmade    = secs;
                
                % send response trigger -- subject accepted an option
                if EEG == 1 && responseTrigNotSent==1
                    sp.sendTrigger(trigger11);
                    responseTrigNotSent=0;
                end
                
            elseif keycode(1,code2) %  
                resp_input  = code2;
                rt          = secs - object_onset;
                answer      = 2; % subject sampled again
                respmade    = secs;
                
                % send response trigger -- subject accepted an option
                if EEG == 1 && responseTrigNotSent==1
                    sp.sendTrigger(trigger12);
                    responseTrigNotSent=0;
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
       
        % BRING FIXATION BACK ON
        Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
        fixation_onset      = Screen('Flip', window, object_offset - slack);    % fixation on, prepare for next trial     
        
        object_offset   = fixation_onset + fixduration +  isi + randperm(jitter*1000,1)/1000 - ifi;  % add jitter here?

        % IF SUBJECT ACCEPTED A FACE/DATE SHOW FFEDBACK AND BREAK SEQUENCE
        % IF SUBJECT CHOSE TO SAMPLE AGAIN, SHOW FIXATION AND MOVE TO THE
        % NEXT FACE
        if answer == 1
            
            % DISPLAY CHOSEN CONTRACT
            feedback_window = Screen('OpenOffscreenWindow', window, windrect);
            Screen('TextSize', feedback_window, textsize);
            Screen('FillRect', feedback_window, grey ,windrect);
            DrawFormattedText(feedback_window, 'Congratulations! This is your new date.', 'center', ycenter-300, white);
            
            % DISPLAY THE CURRENT FACE
            Screen('DrawTexture', feedback_window, textures{thisitem}, [], destrect);     % display the chosen item
            Screen('CopyWindow',feedback_window, window, windrect, windrect);
            object_onset        = Screen('Flip', window, object_offset - slack);    % flip window
            
            % send confidence screen trigger
            if EEG == 1 
                sp.sendTrigger(trigger13)
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
                
                feedback_sampling = Screen('OpenOffscreenWindow',window);
                Screen('TextSize', feedback_sampling, textsize);
                Screen('FillRect', feedback_sampling, grey ,windrect);
                DrawFormattedText(feedback_sampling, 'Oh No! :(', 'center', ycenter-50, white);
                DrawFormattedText(feedback_sampling, 'You you are not allowed to sample again', 'center', ycenter, white);
                Screen('CopyWindow',feedback_sampling, window, windrect, windrect);
                object_onset        = Screen('Flip', window, object_offset - slack);    % flip window
            
                % DISPLAY OFFSET 
                object_offset   = object_onset + feedbacktime - ifi;
                
                % DISPLAY CHOSEN CONTRACT
                feedback_window = Screen('OpenOffscreenWindow', window, windrect);
                Screen('TextSize', feedback_window, textsize);
                Screen('FillRect', feedback_window, grey ,windrect);
                DrawFormattedText(feedback_window, 'Congratulations! This is your new date.', 'center', ycenter-300, white);

                % DISPLAY THE CURRENT FACE
                Screen('DrawTexture', feedback_window, textures{thisitem}, [], destrect);     % display the chosen item
                Screen('CopyWindow',feedback_window, window, windrect, windrect);
                object_onset        = Screen('Flip', window, object_offset - slack);    % flip window

                % send confidence screen trigger
                if EEG == 1 
                    sp.sendTrigger(trigger13)
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
    
    if EEG == 1 
        sp.sendTrigger(trigger101)
    end
    
    chosenitem               = thisitem;

    % save the current trial info
    blocktrials.session      = thisession;
    blocktrials.block        = thisblock;
    blocktrials.trialnumber  = thistrial;
    blocktrials.trialonset   = trialstart;
    blocktrials.sequence     = sequence';
    blocktrials.numsamples   = s;
    blocktrials.chosenitem   = chosenitem;

    set.blocktrials          = blocktrials; % (save the info of this trial) this will go to the main script
    
    WaitSecs(1); % wait two sec before flipping to the next block/

    logs.trials              = trials; % save the samples 

    sublogs                  = fullfile(resfolder,sprintf(logs.trialog,sub,taskname,thisblock,thisession,phase));
    save(sublogs,'logs');

end % end of phase statement 


end