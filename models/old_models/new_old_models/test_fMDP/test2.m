%% define initial variables

% number of trials and draws
Ntrials         = 4;
maxDraws        = 10;

Cc              = 1;
Cw              = 0;
Cs              = -0.025;
q               = 0.6;

% make sequences 
s1              = [1 1 2 2 1 2 1 1 2 1];
s2              = [2 1 1 2 2 1 1 1 2 1];
s3              = [1 1 1 2 2 1 2 1 1 2];
s4              = [1 2 2 1 1 2 2 1 1 1];
drawsequences   = [s1; s2; s3; s4]; % merge sequences

%% run backward induction (calculates utilities and action values)

K               = 3;
reward          = zeros(Ntrials, 1); % initiate reward vector

% loop over trials
for trial = 1:Ntrials
    
    Qsad = zeros(maxDraws, 3); % is that where action values for the 3 options go?
    
    drawSequence = drawsequences(trial, :);
    
    % loop over draws
    for draw = 1:maxDraws
        
        % init utilities
        utility     = zeros(maxDraws, maxDraws+1);
        
        % green bead?
        ng          = sum(drawSequence(1:draw));
        
        % loop over draws backwards
        for drawi       = maxDraws : -1 : (draw + 1) % 10,9,8,7....
            utility_t   = zeros(maxDraws, maxDraws+1);
            futureDraws = drawi - draw;
            ndf         = drawi; % draws 
            
            for greenDraws = 0 : futureDraws
                ngf = ng + greenDraws; % green 
                
                pg = 1/(1 + (q/(1-q))^(ndf-2*ngf));
                pb = 1 - pg;
                
                QG = Cc * pg + Cw * pb;
                QB = Cc * pb + Cw * pg;
                
                if drawi < maxDraws
                    
                    QD = Cs + pb * ((1-q) * utility(ndf+1, ngf+1+1) + (q)*(utility(ndf+1, ngf+1))) +...
                        pg*((q)*utility(ndf+1, ngf+1+1) + (1-q)*(utility(ndf+1, ngf+1)))
                else
                    QD = 0;
                end
                
                
            end % end of greendraws loop
        end % end of backward draw loop 
    end % end of draws loop
end % end of for loop


