
function analyze_model_results_v4;

%assumes my version of spm8 (7/3/2009)
%so far I know that spm_select can't see *.mat files without it

disp('running model');
clear all;
cd E:\beads\formal_subjects\new_models\sub19_cwonly; 
num_subs = 1;
num_conds = 4;
num_trials = num_conds*4;
subject_data = struct('all_sessions_block_types',[],'all_sessions_behavior',[],'all_sessions_sequences',[]);
p1 = 0.8; p2 = 0.6; diff1 = -20; diff2 = -10; amount_sample = -0.25;

    for subject = 1:num_subs; session_data{subject} = spm_select(Inf, 'mesh', sprintf( 'sub %d data', subject')); end;

    for subject = 1:num_subs;  
    
        %clean up after last subject iteration
        clear all_sessions_behavior all_sessions_block_types all_sessions_sequences seq_matrix draws_matrix fitted_output sorted_sess ...
            sorted_trl sorted_types sorted_urns sorted_behavior sorted_subject_data;
        
        %initialize stuff
        onsets = [];
        sess = [];
        trls = [];
        draw_num = [];
        responses = [];
        iterator = 1;
        
        %collapse together session data
        for session = 1:size( session_data{subject} );
    
            clear sequence block_type behavior blue_urn;
            load(session_data{subject}(session,:),'sequence','blue_urn','block_type','behavior');
            
            for trial = 1:size( sequence, 2 );
                
                %if its a green urn sequence then swap the codes so majority is now green
                if blue_urn(trial) == 0;
                    seq_ones = find( sequence{trial} == 1);
                    seq_twos = find( sequence{trial} == 2);
                    sequence{trial}(seq_ones) = 2;
                    sequence{trial}(seq_twos) = 1;
                end;    %ends test for green urn sequences
                
                
                %data which is indexed by sequence/trial
                all_sessions_block_types( iterator ) = block_type( trial );
                all_sessions_sequences{ iterator } = sequence{ trial };
                all_sessions_behavior{ iterator } = behavior{trial};
                iterator = iterator + 1;             
            end;           %ends trials
            
          end;  %ends sessions
            
          subject_data(subject).all_sessions_block_types = all_sessions_block_types;
          subject_data(subject).all_sessions_behavior = all_sessions_behavior;
          subject_data(subject).all_sessions_sequences = all_sessions_sequences;
    
            
            sorted_behavior = [];
            sorted_sess = [];
            sorted_trl = [];
            sorted_types = [];
            sorted_urns = [];
   
    
    %classify cogent output by block_type
    for types = 1:num_conds;
    %for types = 1:1;
        
        %what's this condition?
        switch types;
            case 1; p = p1; cost_diff = diff1; case 2; p = p2; cost_diff = diff1; case 3; p = p1; cost_diff = diff2; case 4; p = p2; cost_diff = diff2;
        end;    %ends switch/case so set block_type parameters p and amounts         
        disp( sprintf('subject %d blocktype: %2.2f p: %2.2f cost_diff: %2.2f sample: %2.2f', subject, types, p,cost_diff, amount_sample)); 
        
        %extract data for this condition
        clear temp type_indices block_batches behavior_batches session_batch trl_batch block_type_batch urn_batch;
        type_indices = find(all_sessions_block_types == types); 
        temp=[];
        for example = 1:size( type_indices, 2);
            
            block_batches(example,:) = all_sessions_sequences{ type_indices(example)  };        %all sessions equences
            behavior_batches{example,:} = all_sessions_behavior{ type_indices(example)}(:,1:3);    %all_sessions_behavior
            temp = [temp; behavior_batches{example,:}];
     
        end;
        
        sorted_behavior = [sorted_behavior; temp];
              
        %%%%%estimate

        [mparams, lla, aQvec] = bayesbeads( block_batches, behavior_batches, subject, types, p, cost_diff);
        fitted_output( types ).params = mparams; 
        fitted_output( types ).lla = lla; 
        fitted_output( types ).aQvec = aQvec;
        
        ftxt = sprintf('fitted_output_%d_%d.mat', subject,types);
        save(ftxt, 'fitted_output', 'behavior_batches', 'block_batches', 'sorted_behavior');
         %figure; %so the subplot will appear in different figure on each run            
%         for i=1:size( fitted_output(types).aQvec, 2 );
%             subplot(6,4,i);
%             h = plot(1:size( fitted_output(types).aQvec{i},1 ), fitted_output(types).aQvec{i});
%             legend('B', 'G', 'D');
%             set(gca, 'Fontname', 'Ariel', 'FontSize', 6);
%             set(h, 'MarkerSize',6, 'marker', 'o');
%         end;
%         title( sprintf('subject %d condition %d',subject, types) );
        
    end;    %ends loop through block types
  
  end;  %ends loop through subjects
  
  disp('audi5000');
    
%%    
function [mparams, lla, aQvec] = bayesbeads(sequence, choiceVec, subject, types, p, cost_diff)

G = 1;
B = 2;
S = 3;

Ntrials  = 2;
maxDraws = 10;

if subject == -1
    %%% sets up outcomes for each bead draw
    sequence = (rand(Ntrials, maxDraws) > 0.5) + 1;
    choiceVec = setChoiceVec(sequence); 
end
    
alpha = 1; %%% softmax function which converts values to probabilities
%q     = 0.8; %% probability ratio of green vs. blue beads
q=p;
Cw    = -50; %% cost of being wrong
%Cw = cost_diff;
Cs    = -0.25;  %% cost of a sample

%%% compile parameters into vector
params      = [Cw ];

fixedParams = [alpha; q; Cs];

%%% switch which controls what estimatelikelihood does
findPick = 0;

% estimateLikelihood(params, sequence, choiceVec, q, findPick);

%%% determine optimal stopping point across a range of parameter values
% spaceDat = explorePickSpace();
% 
% save('spaceDat.mat', 'spaceDat');

%%% optimizing parameters for invidual subjects

options = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);

initialJitter = [ 5; -5; 5; -5];
%                  [ 5   0.10;
%                  -5   0.10;
%                   5  -0.10;
%                  -5  -0.10];
                
 llaMin = Inf;
%  for startValue = 1 : 4
    
    %startParam = params + initialJitter(startValue, :)'; 
    startParam = params;
    
    [mparams, lla] = fminsearch(@(params) estimateLikelihood(params, sequence, choiceVec, fixedParams, findPick),startParam, options);
    
    if lla < llaMin
        llaMin = lla;
        minParams = mparams;
    end
    
 %end

fprintf('ll %.3f\n', lla);

% ftxt = sprintf('subjectParams_%d_%d.mat', subject,types);
% save(ftxt, 'minParams', 'llaMin');

[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, sequence, choiceVec, fixedParams, findPick);

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this function looks at optimal stopping points for different parameter
%%% values, q, Cw and Cs
function spaceDat = explorePickSpace()

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
            
            params = [Cw, Cs, alpha];

            %%% run for 50 different sets of sequences, because sequences
            %%% are stochastic themselves
            for boot_sequence = 1 : 50

                %%% generate a set of random sequences, given q
                sequence = (rand(Ntrials, maxDraws) > q) + 1;

                %%% setup choicevec so analysis runs to end of draws
                for set = 1 : Ntrials
                %     draw = ceil((maxDraws-4)*rand(1,1));
                    draw = maxDraws;
                    choiceVec{set} = zeros(draw, 3);
                    choiceVec{set}((1:(end-1)),3) = 1;
                    if mean(sequence(set, 1:draw)) > 1.5
                        choiceVec{set}(end, 2) = 1;
                    else
                        choiceVec{set}(end, 1) = 1;
                    end
                end
                
                %%% determine where in sequence of draws optimal subject
                %%% would stop
                [ll, pickTrial] = estimateLikelihoodf(params, sequence, choiceVec, fixedParams, findPick);
                
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ll] = estimateLikelihood(params, sequence, setData, fixedParams, findPick)

[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setData, fixedParams, findPick);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setData, fixedParams, findPick)

%%% extract parameters from parameter vector
Cw = params(1);
%Cs = params(2);
%q  = params(2);

%Cs    = fixedParams(1);
alpha = fixedParams(1);
q = fixedParams(2);
Cs = fixedParams(3);

%%% length of sequence
lseq    = size(sequence, 2);

%%% number of sequences
nblocks = size(setData, 1);

%%% intialize log likelihood to zero
ll = 0;

for block = 1 : nblocks
    
    %%% choices of subject for this sequence of draws
    choiceVec = setData{block};

    %%% number of choices for this sequence
    nchoices = size(choiceVec, 1);
    
    %%% initially nd (draws) == 0 and ng (green marbles) == 0
    ng = 0;
    nd = 0;
    
    %%% Qvec is values of each action
    dQvec = [];
    
    %%% corresponding probabilities generated with softmax and alpha
    ddec  = []; 
        
%     if isempty(find( mean( choiceVec') == 0 )) == 0; 
%         disp(sprintf('missing response skipping sequence %d for this type', block)); continue; 
%     end;

    %%% loop over draws for this sequence of draws
    for draw = 1 : nchoices

        %%% check if we got a green or blue marble
        if sequence(block, draw) == 1
            ng = ng + 1;
        end

        %%% alwasy increment draws
        nd = nd + 1;

        %%% compute values of each action
        [v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cs);
        
        %%% keep track of values across sequnce of draws
        dQvec(draw, 1:length(Qvec)) = Qvec;
        
        %%% keep track of choice probabilities across sequence 
        ddec (draw, 1:length(d))    = d;
        
        %%% Nick, add this line to your code and also return it.
        aQvec{block}(draw, 1:length(Qvec)) = Qvec;
        
        %%% if trying to determine optimal stopping position (findpick ==
        %%% 1) then see if we should stop
        if findPick == 1 & draw < nchoices & (Qvec(1) > Qvec(3) | Qvec(2) > Qvec(3))
            pickTrial(block) = draw;
            break
        elseif findPick == 1 & draw == nchoices
            pickTrial(block) = draw;
        end

        %%% if at end of sequence d is two element and have to tack on for
        %%% syntax purposes
        if draw == lseq
            d = [d; 0];
        end
        
        if choiceVec(draw,:)*d == 0;
            fprintf('missing data sequence %d   ', block);
            choiceVec(draw,:) = [1/3 1/3 1/3];
        end
            

        %%% update log likelihood
        try
            ll = ll - log(choiceVec(draw,:)*d);
        catch
            fprintf('');
        end
        
        %fprintf( 'draw: %d ll: %.2f', nd, ll);
        
%         if findPick == 0
% 
%             if length(Qvec) == 3
%                 fprintf('%d\t', block);
%                 for di = 1 : 3
%                     fprintf('%.2f/%.2f\t', Qvec(di), d(di));
%                 end
%                 fprintf('%d pg %.2f\n', sequence(block, draw), PG(q, nd, ng));
%             end
%              
%          end

    end
    
%     subplot(2,2,block);
%     plot(dQvec);
%     legend('G', 'B', 'S');

%     subplot(6,4,block);
%     h = plot(1:size(dQvec,1), dQvec);
%     legend('B', 'G', 'D');
%     set(gca, 'Fontname', 'Ariel', 'FontSize', 6);
%     set(h, 'MarkerSize',6, 'marker', 'o');
%     
end

if findPick == 0
    pickTrial = [];
    fprintf('ll %.2f Cw %.1f p %.2f Cs %.2f alpha %.2f\n', ll, Cw, q, Cs, alpha);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v] = vVal(q, nd, ng, alpha, lseq, Cw, Cs)

[v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cs)

%%% computes probability that we are drawing from green urn
pg = PG(q, nd, ng);
%%% probability that we are drawing from predominantly blue urn
pb = 1 - pg;

%%% cost of choosing Green, at this point
QG = Cw*pb;
%%% cost of choosing blue at this point
QB = Cw*pg;


if nd + 1 <= lseq
    try
        
        %%% compute value of next state given that we draw a green
        val11 = vVal(q, nd+1, ng+1, alpha, lseq, Cw, Cs);
        %%% compute value of next state given that we draw a blue
        val10 = vVal(q, nd+1, ng, alpha, lseq, Cw, Cs);
        
        %%% redundant with above
%         val00 = vVal(q, nd+1, ng+1, alpha, lseq, Cw, Cs);
%         val01 = vVal(q, nd+1, ng, alpha, lseq, Cw, Cs);

        %%% Value of action is cost to sample plus expected value of next
        %%% state
        QS = Cs + pg*(val11*q     + val10*(1-q)) +...
                  pb*(val11*(1-q) + val10*q);
              
        %%% compile action values into vector
        Qvec = [QG; QB; QS]; 
        if nd == 1
           fprintf(''); 
        end
    catch
        fprintf('');
    end
else
    Qvec = [QG; QB];
end

try sum(exp(alpha*Qvec));
catch; fprintf('');
end

if sum(exp(alpha*Qvec)) == 0
    fprintf('');
end

%%% softmax function to convert values to action probabilities
d = exp(alpha*Qvec)./sum(exp(alpha*Qvec));

%%% average value of this state if we take actions with probability d
v = d'*Qvec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function choiceVec = setChoiceVec(sequence)

%%% Ntrials is total number of bead draw sequences
for set = 1 : Ntrials
%     draw = ceil((maxDraws-4)*rand(1,1));
    draw = maxDraws;
    %%% initialize matrix
    choiceVec{set} = zeros(draw, 3);
    %%% set choice to draw again up to last bead
    choiceVec{set}((1:(end-1)),3) = 1;
    %%% solve task rationally at last bead draw of sequence
    if mean(sequence(set, 1:draw)) > 1.5
        choiceVec{set}(end, 2) = 1;
    else
        choiceVec{set}(end, 1) = 1;
    end
end

