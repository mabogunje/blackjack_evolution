%{
    author: Damola Mabogunje
    contact: damola@mabogunje.net
    summary: HW6 Blackjack Player Evolution
             An evolutionary approach to learning to win a game of blackjack.
%}

%% Initializing Variables
SELECTION_CRITERIA = 0.50; % Min Fitness Selection critieria
MUTATION_FREQUENCY = 0.10; % Default frequency of mutation
TERMINAL_CRITERIA = 3.00; % Required improvement over Dealer
STATS = []; % For recordings key statistics

% Game Data
GOAL = 21; % Blackjack!
GAMES = 100; % # of games to base fitness on
STATES = [4:12 12.5:0.5:20.5 21 22]; % Floats represent soft states, 22 is BUST
NUM_STATES = numel(STATES);
sort(STATES);

% Player Policies (Player [P], Dealer [D])
Ppolicy = randi([0,1], NUM_STATES, 1);
Dpolicy = randi([0,1], NUM_STATES, 1);

% Dealer stands below 17 (trivial)
%Dpolicy = zeros(NUM_STATES, 1);
%Dpolicy(19:end) = 1;


% Force Stand on BUST (Prevents infinite loop)
Ppolicy(end) = 0;
Dpolicy(end) = 0;

% Start States
Pstate = STATES(randperm(NUM_STATES, 1));
Dstate = STATES(randperm(NUM_STATES, 1));

% Player Records
best = 0.0001; % ::KLUDGE:: Avoids div by zero in later computations
Precord = zeros(GAMES, NUM_STATES);
% Drecord is unassigned for now

%% Begin
generations = 1;

fprintf('Goal: %d%% Improvement of Player policy over Dealer policy\n', TERMINAL_CRITERIA*100);
disp('Dealer, Player');
fprintf('%d,\t%d \n\n', Dpolicy, Ppolicy);
disp('Evolving...');
while(best < TERMINAL_CRITERIA)
    
    % Simulate Games
    for simulation = 1 : GAMES
        for state  = 1 : NUM_STATES
            Dstate = STATES(state);
            Pstate = STATES(state);
        
            %disp([Dstate, Dpolicy(state), Ppolicy(state)]);

            while(Dpolicy(ceil(Dstate-STATES(1)+1)) > 0) % Dealer hit
                transitions = STATES(find(STATES==Dstate)+1:NUM_STATES);
                    
                if(~isempty(transitions))
                    Dstate = transitions(randperm(numel(transitions), 1));
                else
                    break;
                end
            end

            while(Ppolicy(ceil(Pstate-STATES(1)+1)) > 0) % Player hit
                transitions = STATES(find(STATES==Pstate)+1:NUM_STATES);

                if(~isempty(transitions))
                    Pstate = transitions(randperm(numel(transitions), 1));
                else
                    break;                
                end
            end                
            
            %disp([Dstate, Pstate]);
        
            % Record Policy Outcomes
            if(Pstate > GOAL || (Dstate < GOAL && Pstate < Dstate))
                Precord(simulation, state) = 0;
            else
                Precord(simulation, state) = 1;
            end
        end    
    end

    % Drecord can now be assigned as game is adversarial
    Drecord = ~Precord;

    %% Fitness Evaluation
    Pfit = sum(Precord)/GAMES;
    Dfit = sum(Drecord)/GAMES;
    Psum = cumsum(sum(Precord));
    Dsum = cumsum(sum(Drecord));
    newBest = ((Psum(end) - Dsum(end)) / Dsum(end)) - 1;

    %% Inheritance
    Pchild = Ppolicy;
    Dchild = Dpolicy;

    %% Selection
    Pgenes = Pfit - Dfit;
    Dgenes = Dfit - Pfit;

    %% Crossover
    % We psuedo-coevolve the Dpolicy via Dchild 
    % for the benefit of improved crossover
    
    fitGenes = find(Pgenes >= SELECTION_CRITERIA);
    Dchild(fitGenes) = Ppolicy(fitGenes);

    fitGenes = find(Dgenes >= SELECTION_CRITERIA);
    Pchild(fitGenes) = Dchild(fitGenes);
    
    %% Mutation
    if((newBest - best) > 0)
        mutRate = MUTATION_FREQUENCY;
    else % Adjust according to how much we sucked!
        mutRate = abs(newBest - best);
        
        % If we reaally sucked, remember to max at 100% mutation
        if(mutRate > 1)
            mutRate = 1;
        end
    end

    nummut = ceil(mutRate * (NUM_STATES-1)); % -1 excludes BUST which always stands
    mutants = randi(NUM_STATES, nummut, 1);
    Pchild(mutants) = ~Pchild(mutants);

    %% Evolution 
    % Evolve Player if we see improvement
    if(newBest > best)
        best = newBest;
        Ppolicy = Pchild;
    else % Otherwise mutate current policy
        Ppolicy(mutants) = ~Ppolicy(mutants);
    end
    
    generations = generations + 1;
    
    %% Results
    crossedover = numel(fitGenes);
    mutations = numel(mutants); % Alias for understanding results
    geneFitness = Psum(end) / (GAMES*generations); % Also average fitness
    
    results = [generations, geneFitness, newBest*100, crossedover, mutations, TERMINAL_CRITERIA*100];
    STATS = [STATS; results];
   
end

%% Output
headers = {'State' 'Dealer', 'Player'};
columns = {'Generation', 'Avg. Gene Fitness', 'Policy Fitness', 'Crossovers', 'Mutations', 'Goal'};
data = '\t%i,\t\t %+.2f,\t\t\t %+.2f,\t\t %i,\t\t %i,\t\t %+.2f \n';

disp('Done!');

disp('Statistics');
disp(columns)
fprintf(data, STATS');
disp(sprintf('\n\n'));
disp('Dealer, Player');
fprintf('%i,\t%i \n', Dpolicy, Ppolicy);

%% Export to CSV
EXPORT_FILE = 'BPE_complex.csv';
csvwrite_with_headers(EXPORT_FILE, STATS, columns);



