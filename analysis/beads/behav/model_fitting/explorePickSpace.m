function spaceDat = explorePickSpace()

% this function looks at optimal stopping points for different parameter values, q, Cw and Cs

% define parameters
ntrials     = 1; 
maxdraws    = 10;
boots       = 5;
alpha       = 1;
findpick    = 1;
qctr        = 1; % what is this?

for q = 0.6 : 0.2 : 0.8
    
    wctr = 1; % what is this?
    
    for costwrong = -100 : 20 :0
        
        fprintf('q %.1f cost wrong %d\n', q, costwrong);
        
        sctr = 1;
        
        for costsample = -1 : 1 : 0
            
            params = [costwrong, costsample, alpha];
            
            % run bootstraps on sets of random sequences 
            for bootsequence = 1:boots
                
                % create a random sequence (with either 0.6 or 0.8
                % probability)
                random_seq = (rand(ntrials, maxdraws) > q) + 1;
                
                % set up choice vector (what is that?)
                for set = 1: ntrials
                    
                    draw = maxdraws;  % number of draws
                    choicevec{set} = zeros(draw, 3); % set the matrix 
                    choicevec{set}((1:(end-1)),3) = 1; % choices for drawing again up to last bead 
                    
                    if mean(random_seq(set, 1:draw)) > 1.5
                        choicevec{set}(end, 2) = 1
                        
                    else
                        choicevec{set}(end, 1) = 1
                    end
                    
                end % end of choice vector loop
                
                % determine where in sequence of draws optimal subject would stop
                [ll, picktrial] = estimateLikelihoodf(params, random_seq, choicevec, fixedparams, findpick);
                
                %%% keep track of point in sequence where subject stopped
                bootpick(bootsequence) = picktrial;
                
            end % end of bootstrap loop
            
            %%% how do parameters affect optimal stopping point
            spacedat(qctr, wctr, sctr) = mean(bootpick);
            
            sctr = sctr + 1;
            
        end % end of cost sample for loop
        
        wctr = wctr + 1;
        
    end % end of cost wrong for loop
    
    qctr = qctr + 1;

end % end of q for loop

return