% Function to simulate the participant's responses using the learned policy
function [nDraws, urnChoice] = simulateParticipant(beads, policy, blueRatio)
    nBeads = length(beads);
    state = 1;
    nDraws = 0;

    % Loop until the participant makes an urn choice
    while true
        % Compute the current urn color
        urnColor = 2 - state;

        % Choose the action according to the optimal policy
        action = policy(state);

        % Update the state and the number of draws for the chosen action
        if action == 1
            % Guess the urn color
            guess = 2 - urnColor;
            isCorrect = guess == beads(nDraws+1);
            if isCorrect
                urnChoice = urnColor;
            else
                urnChoice = 2 - urnColor;
            end
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
end
