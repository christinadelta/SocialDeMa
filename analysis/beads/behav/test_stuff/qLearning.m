% Function to learn the optimal policy using Q-learning
function [Q, policy,nDraws] = qLearning(beads, blueRatio, rewardCorrect, rewardIncorrect, costSample)

nBeads = length(beads);
Q = zeros(2, 2);
alpha = 0.2;
gamma = 0.2;
epsilon = 0.8;
delta = inf;
iter = 1;

while delta > 1e-6 && iter <= nBeads
    % Initialize the state and the number of draws
    state = 1;
    nDraws = 0;
    Qprev = Q;

    % Loop until the participant makes an urn choice
    while true
        % Choose the action according to the epsilon-greedy policy
        if rand < epsilon
            action = randi(2);
        else
            [~, action] = max(Q(state, :));
        end

        % Update the state and the number of draws for the chosen action
        if action == 1
            % Guess the urn color
            guess = 2 - state;
            isCorrect = guess == beads(iter);
            if isCorrect
                reward = rewardCorrect - (costSample * nDraws);
            else
                reward = rewardIncorrect - (costSample * nDraws);
            end
            Q(state, action) = Q(state, action) + alpha * (reward - Q(state, action));
            iter = iter + 1;
            break;
        else
            % Sample more beads
            nDraws = nDraws + 1;
            if rand < blueRatio
                state = 1;
            else
                state = 2;
            end
        end
    end

    % Compute the change in Q-values from the previous iteration
    delta = max(abs(Q(:) - Qprev(:)));
end

% Compute the optimal policy from the learned Q-values
[~, policy] = max(Q, [], 2);


end
