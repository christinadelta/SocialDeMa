function spaceDat = explorePickSpace()

% this function looks at optimal stopping points for different parameter values, q, Cw and Cs

% define parameters
ntrials     = 1; 
maxdraws    = 10;
alpha       = 1;
findpick    = 1;
qctr        = 1;

for q = 0.6 : 0.2 : 0.8
    
    wctr = 1;
    
    for costworng = -100 : 20 :0
        
        fprintf('q %.1f cost wrong %d\n', q, costworng);
        
        sctr = 1;
        
        for costsample = -1 : 1 : 0
            
            params = [costwrong, costsample, alpha];
            
            
        end % end of cost sample for loop
        
    end % end of cost wrong for loop

end % end of q for loop

















return