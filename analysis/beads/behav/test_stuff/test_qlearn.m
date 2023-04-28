

clear all
clc

%%
% Set the fixed parameters
nSequences = 10;
nBeads = 10;
blueRatio = 0.6;
rewardCorrect = 10;
rewardIncorrect = -10;
costSample = 0.25;

%%
% Initialize the arrays to store the results
drawCount = zeros(nSequences, 1);
urnChoices = zeros(nSequences, 1);

%%
% Loop over the sequences
for i = 1:nSequences
    % Generate the sequence of beads
    beads = generateBeadSequence(nBeads, blueRatio);
    
    % Run the Q-learning algorithm to learn the optimal policy
    [Q, policy,nDraws(i)] = qLearning(beads, blueRatio, rewardCorrect, rewardIncorrect, costSample);
    
    % Simulate the participant's responses using the learned policy
    [drawCount(i), urnChoices(i)] = simulateParticipant(beads, policy, blueRatio);
end

%%
% Display the results
disp("Number of draws:");
disp(nDraws);
disp("Urn choices:");
disp(urnChoices);

