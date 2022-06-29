%%% testing synthetic data and bootstarping 

Ntrials  = 1;
maxDraws = 10;
    
alpha = 1;

findPick = 1;

qCtr = 1;
for q = 0.6 : 0.2 : 0.8
        
    wCtr = 1;
    for Cw = -100 : 20 : 0
        
        fprintf('q %.1f Cw %d\n', q, Cw);
        
        sCtr = 1;
        for Cs = -6 : 1 : 0
            
            fixedParams = [alpha; q; Cs];
            params = Cw;
            % params = [Cw, Cs, alpha];

            %%% run for 50 different sets of sequences, because sequences
            %%% are stochastic themselves
            for boot_sequence = 1 : 20

                %%% generate a set of random sequences, given q
                sequence = (rand(Ntrials, maxDraws) > q) + 1;

                %%% setup choicevec so analysis runs to end of draws
                for set = 1 : Ntrials
                %     draw = ceil((maxDraws-4)*rand(1,1));
                    draw = maxDraws; % number of draws
                    choicevec{set} = zeros(draw, 3); % set the matrix 
                    choicevec{set}((1:(end-1)),3) = 1; % choices for drawing again up to last bead 
                    
                    if mean(sequence(set, 1:draw)) > 1.5
                        choicevec{set}(end, 2) = 1;
                    else
                        choicevec{set}(end, 1) = 1;
                    end
                end
                
                %%% determine where in sequence of draws optimal subject
                %%% would stop
                [ll, pickTrial] = estimateLikelihoodf(params, sequence, choicevec, fixedParams, findPick);
                
                %%% keep track of point in sequence where subject stopped
                bootPick(boot_sequence) = pickTrial;

            end
            
            %%% how do parameters affect optimal stopping point
            spaceDat(qCtr, wCtr, sCtr) = mean(bootPick);
            
            sCtr = sCtr + 1;
            
        end
        wCtr = wCtr + 1;
    end
    qCtr = qCtr + 1;
end