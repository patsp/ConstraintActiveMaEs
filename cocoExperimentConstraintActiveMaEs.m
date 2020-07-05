% This script runs the Constraint Active-MA-ES for BUDGET_MULTIPLIER*DIM
% function evaluations on the 'bbob-constrained' suite.
%
% It is based on the exampleexperiment.m provided by the BBOB COCO
% framework.

more off; % to get immediate output in Octave

name = 'ConstraintActiveMaEs';
firstFunction = 1;
lastFunction = 48;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%
BUDGET_MULTIPLIER = 100000; % algorithm runs for BUDGET_MULTIPLIER*dimension funevals
NUM_OF_INDEPENDENT_RESTARTS = 0; % max. number of independent algorithm
% restarts; if >0, make sure that the
% algorithm is not always doing the same thing
% in each run (which is typically trivial for
% randomized algorithms)

%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Experiment    %
%%%%%%%%%%%%%%%%%%%%%%%%%
suite_name = 'bbob-constrained'; % works for 'bbob' as well
%suite_name = 'bbob';
observer_name = 'bbob';
observer_options = strcat('result_folder: ', name, '_on_', ...
                          suite_name, '_f', num2str(firstFunction), ...
                          '_', num2str(lastFunction), ...
                          ' algorithm_name: ', name, ...
                          ' algorithm_info: ', name, ...
                          ' target_precision: 1e-8');

% dimension 40 is optional:
suite = cocoSuite(suite_name, 'instances: 1-15', 'dimensions: 2,3,5,10,20,40');
observer = cocoObserver(observer_name, observer_options);

% set log level depending on how much output you want to see, e.g. 'warning'
% for fewer output than 'info'.
cocoSetLogLevel('info');

% keep track of problem dimension and #funevals to print timing information:
printeddim = 1;
doneEvalsAfter = 0; % summed function evaluations for a single problem
doneEvalsTotal = 0; % summed function evaluations per dimension
printstring = '\n'; % store strings to be printed until experiment is finished

nInstances = 15;
nTargetProblemBegin = firstFunction;
nTargetProblemEnd = lastFunction;
prevDimension = 0;

cnt = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Experiment        %
%%%%%%%%%%%%%%%%%%%%%%%%%
while true

    % get next problem and dimension from the chosen suite:
    problem = cocoSuiteGetNextProblem(suite, observer);
    if ~cocoProblemIsValid(problem)
        break;
    end
    dimension = cocoProblemGetDimension(problem);

    if dimension ~= prevDimension
      cnt = 1;
    end
    prevDimension = dimension;

    good = false;
    if (nTargetProblemBegin - 1) * nInstances < cnt && ...
       cnt < nTargetProblemEnd * nInstances + 1
        good = true;
    end
    if ~good
      cnt = cnt + 1;
      continue;
    end

    % printing
    if printeddim < dimension
      if printeddim > 1
        elapsedtime = toc;
        printstring = strcat(printstring, ...
            sprintf('   COCO TIMING: dimension %d finished in %e seconds/evaluation\n', ...
            printeddim, elapsedtime/double(doneEvalsTotal)));
        tic;
      end
      doneEvalsTotal = 0;
      printeddim = dimension;
      tic;
    end

    % restart functionality: do at most NUM_OF_INDEPENDENT_RESTARTS+1
    % independent runs until budget is used:
    i = -1; % count number of independent restarts
    while (BUDGET_MULTIPLIER*dimension > (cocoProblemGetEvaluations(problem) + ...
                                          cocoProblemGetEvaluationsConstraints(problem)))
        i = i+1;
        if (i > 0)
            fprintf('INFO: algorithm restarted\n');
        end
        doneEvalsBefore = cocoProblemGetEvaluations(problem) + ...
                          cocoProblemGetEvaluationsConstraints(problem);

        % start algorithm with remaining number of function evaluations:
        runConstraintActiveMaEsOnCocoProblem(...
                            problem,...
                            cocoProblemGetSmallestValuesOfInterest(problem),...
                            cocoProblemGetLargestValuesOfInterest(problem), ...
                            BUDGET_MULTIPLIER*dimension - doneEvalsBefore);

        % check whether things went wrong or whether experiment is over:
        doneEvalsAfter = cocoProblemGetEvaluations(problem) + ...
                         cocoProblemGetEvaluationsConstraints(problem);

        if cocoProblemFinalTargetHit(problem) == 1 ||...
             doneEvalsAfter >= BUDGET_MULTIPLIER * dimension
            break;
        end
        if (doneEvalsAfter == doneEvalsBefore)
            fprintf('WARNING: Budget has not been exhausted (%d/%d evaluations done)!\n', ...
                doneEvalsBefore, BUDGET_MULTIPLIER * dimension);
            break;
        end
        if (doneEvalsAfter < doneEvalsBefore)
            fprintf('ERROR: Something weird happened here which should not happen: f-evaluations decreased');
        end
        if (i >= NUM_OF_INDEPENDENT_RESTARTS)
            break;
        end
    end

    doneEvalsTotal = doneEvalsTotal + doneEvalsAfter;

    cnt = cnt + 1;
end

elapsedtime = toc;
printstring = strcat(printstring, ...
    sprintf('   COCO TIMING: dimension %d finished in %e seconds/evaluation\n', ...
    printeddim, elapsedtime/double(doneEvalsTotal)));
fprintf(printstring);

cocoObserverFree(observer);
cocoSuiteFree(suite);
