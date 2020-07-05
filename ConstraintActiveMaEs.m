% Copyright 2018 The Constraint Active-MA-ES Authors. All Rights Reserved.
%
% This file is part of Constraint Active-MA-ES.
%
% Constraint Active-MA-ES is free software: you can redistribute it
% and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation, either version 3 of
% the License, or (at your option) any later version.
%
% Constraint Active-MA-ES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Constraint Active-MA-ES. If not, see
% <http://www.gnu.org/licenses/>.

function [x, info] = ConstraintActiveMaEs(objectiveFun, constraintFun, config)
  % [x, info] = ConstraintActiveMaEs(objectiveFun, constraintFun, config)
  %    Attempts to minimize objectiveFun(x)
  %    subject to y(1) <= 0, ..., y(nConstraints) <= 0
  %    where y = constraintFun(x) and y and x are expected to be column vectors.
  %
  %    The parameters are:
  %        objectiveFun is expected to be a handle to a fitness function. It
  %            should take one argument (column vector x), evaluate the
  %            objective function on x and return the result.
  %        constraintFun is expected to be a handle to a constraint function. It
  %            should take one argument (column vector x), evaluate the
  %            constraints on x and return a column vector y that contains x
  %            evaluated on all constraints.
  %        config is expected to be a struct containing configuration parameters
  %            config.dimension: the dimension of the problem.
  %            config.xInit: the value for the initial populations centroid.
  %            config.observer: handle to a callback function (can for example
  %                             be used for logging purposes).
  %            config.maxEvals: the sum of the total number of objective
  %                             function and constraint evaluations that can
  %                             be used.
  %
  %    The return values are:
  %        x is the final population's centroid.
  %        info is an information struct containing:
  %            info.aBest: the best individual found (struct containing
  %                        the best found parameter vector x and its
  %                        objective function value f)
  %            info.terminationCriterionStr: termination criterion description
  %            info.terminationCriterion: termination criterion id
  %            info.numGenerations: number of generations
  %            info.nEvals: number of objective function and constraint
  %                         evaluations
  %
  %    This is the implementation of the Constraint Active-MA-ES from [3].
  %    The idea is to extend the MA-ES described in [1] with the idea
  %    of using the Active Covariance Update for constraint handling
  %    proposed in [2] for the case of the (1+1)-CMA-ES.
  %    [1]
  %    Simplify Your Covariance Matrix Adaptation Evolution Strategy
  %    Beyer, Hans-Georg. and Sendhoff, Bernhard
  %    (2017)
  %    [2]
  %    A (1+1)-CMA-ES for Constrained Optimisation
  %    Arnold, Dirk and Hansen, Nikolaus
  %    (2012)
  %    [3]
  %    A multi-recombinative active matrix adaptation evolution
  %    strategy for constrained optimization
  %    Spettel, Patrick and Beyer, Hans-Georg

  dimension = config.dimension;

  nFitnessEvaluations = 0;
  nConstraintEvaluations = 0;

  x = zeros(dimension, 1);
  if isfield(config, 'xInit')
    x = config.xInit;
  end
  info.aBest.x = x;
  info.aBest.f = objectiveFun(x);
  aBestG = 0;

  nConstraints = size(constraintFun(x), 1);
  nConstraintEvaluations = nConstraintEvaluations + 1;

  epsilon = 1e-20;
  gStop = 100000;
  sigmaStop = 1e-10;
  lambda = 4 + floor(3 * log(dimension));
  mu = floor(lambda / 2);
  weightsTemp = log((lambda + 1) / 2) - log((1:lambda)');
  muEff = (sum(weightsTemp(1:mu))) ^ 2 / ...
           (sum(weightsTemp(1:mu) .^ 2));
  muEffNeg = (sum(weightsTemp((mu+1):lambda))) ^ 2 / ...
               (sum(weightsTemp((mu+1):lambda) .^ 2));
  cSigma = (muEff + 2) / (dimension + muEff + 5);
  cC = (4 + muEff / dimension) / ...
        (dimension + 4 + 2 * muEff / dimension);
  c1 = 2 / ((dimension + 1.3) ^ 2 + muEff);
  cMu = min(1 - c1, ...
             2 * ((muEff - 2 + 1 / muEff) / ...
                  ((dimension + 2) ^ 2 + 2 * muEff / 2)));
  alphaMuNeg = 1 + c1 / cMu;
  alphaMuEffNeg = 1 + (2 * muEffNeg) / (muEff + 2);
  alphaPosdefNeg = (1 - c1 - cMu) / ...
                     (dimension * cMu);
  weights = zeros(lambda, 1);
  positiveWeightIndices = weightsTemp >= 0;
  weights(positiveWeightIndices) = ...
    (1 / sum(weightsTemp(positiveWeightIndices))) * ...
    weightsTemp(positiveWeightIndices);
  negativeWeightIndices = weightsTemp < 0;
  weights(negativeWeightIndices) = ...
    (min([alphaMuNeg, alphaMuEffNeg, alphaPosdefNeg]) / ...
     sum(abs(weightsTemp(negativeWeightIndices)))) * ...
    weightsTemp(negativeWeightIndices);
  assert(sum(weights > 0) == floor(lambda / 2));
  sigma = 1;
  cV = 1 / (dimension + 2);
  M = eye(dimension);
  V = zeros(dimension, nConstraints);
  s = zeros(dimension, 1);
  gLag = 50 * dimension;
  beta = 0.1 / (dimension + 2);

  g = 0;
  while true
    if isfield(config, 'observer')
      config.observer(x, info.aBest.x, ...
                      nConstraintEvaluations + nFitnessEvaluations);
    end

    offspringPop = cell(lambda);
    fitnesses = zeros(lambda, 1);

    violatedConstraints = false(nConstraints, 1);
    nFeasible = 0;
    for k = 1:lambda
      offspring.z = randn(dimension, 1);
      offspring.d = M * offspring.z;
      offspring.x = x + sigma * offspring.d;

      offspring.constraintValues = constraintFun(offspring.x);
      nConstraintEvaluations = nConstraintEvaluations + 1;
      offspring.feasible = all(offspring.constraintValues <= epsilon);
      if offspring.feasible
        nFeasible = nFeasible + 1;
        offspring.fitness = objectiveFun(offspring.x);
        nFitnessEvaluations = nFitnessEvaluations + 1;
        if offspring.fitness < info.aBest.f
          info.aBest.x = offspring.x;
          info.aBest.f = offspring.fitness;
          aBestG = g + 1;
        end
      else
        offspring.fitness = inf;
        for j = 1:nConstraints
          if offspring.constraintValues(j) > epsilon
            V(:, j) = (1 - cV) * V(:, j) + cV * offspring.z;
            violatedConstraints(j) = true;
          end
        end
      end

      fitnesses(k) = offspring.fitness;
      offspringPop{k} = offspring;
    end

    [~, indicesSorted] = sort(fitnesses);

    if nFeasible >= mu
      zW = zeros(dimension, 1);
      dW = zeros(dimension, 1);
      dzW = zeros(dimension, dimension);
      for k = 1:nFeasible
        z = offspringPop{indicesSorted(k)}.z;
        d = offspringPop{indicesSorted(k)}.d;
        if k <= mu
          zW = zW + weights(k) * z;
          dW = dW + weights(k) * d;
        end
        dzW = dzW + weights(k) * d * z';
      end

      xTmp = x + sigma * dW;
      feasible = all(constraintFun(xTmp) <= epsilon);
      nConstraintEvaluations = nConstraintEvaluations + 1;
      if feasible
        x = xTmp;

        f = objectiveFun(x);
        nFitnessEvaluations = nFitnessEvaluations + 1;

        if f < info.aBest.f
          info.aBest.x = x;
          info.aBest.f = f;
          aBestG = g + 1;
        end

        s = (1 - cSigma) * s + ...
            sqrt(muEff * cSigma * (2 - cSigma)) * zW;
        M = (1 - (c1 / 2) - (cMu / 2)) * M + ...
            (c1 / 2) * (M * s) * s' + ...
            (cMu / 2) * dzW;
        sigma = sigma * exp((cSigma / 2) * ...
                            (((s' * s) / dimension) - 1));
      else
        fx = objectiveFun(x);
        nFitnessEvaluations = nFitnessEvaluations + 1;
        bestOffspring = offspringPop{indicesSorted(1)};
        if bestOffspring.fitness < fx
          x = bestOffspring.x;

          f = objectiveFun(x);
          nFitnessEvaluations = nFitnessEvaluations + 1;

          if f < info.aBest.f
            info.aBest.x = x;
            info.aBest.f = f;
            aBestG = g + 1;
          end

          s = (1 - cSigma) * s + ...
              sqrt(muEff * cSigma * (2 - cSigma)) * bestOffspring.z;
          M = (1 - (c1 / 2) - (cMu / 2)) * M + ...
              (c1 / 2) * (M * s) * s' + ...
              (cMu / 2) * (M * bestOffspring.z) * bestOffspring.z';
          sigma = sigma * exp((cSigma / 2) * ...
                              (((s' * s) / dimension) - 1));
        end
      end
    end

    for j = 1:nConstraints
      if violatedConstraints(j)
        M = M - beta * (M * V(:, j)) * V(:, j)';
      end
    end

    g = g + 1;
    if nFitnessEvaluations + nConstraintEvaluations > config.maxEvals || ...
       g > gStop || sigma < sigmaStop || g - aBestG >= gLag
      break;
    end
  end

  info.terminationCriterionStr = '';
  info.terminationCriterion = 0;
  if g > gStop
    info.terminationCriterionStr = ...
      ['Maximum number of generations reached, g > gStop ', ...
       sprintf('(%d > %d)', g, gStop)];
    info.terminationCriterion = 1;
  elseif sigma < sigmaStop
    info.terminationCriterionStr = ...
      ['Mutation step size decreased below threshold,', ...
       ' sigma < sigmaStop ', ...
       sprintf('(%f < %f)', sigma, sigmaStop)];
    info.terminationCriterion = 2;
  elseif g - aBestG >= gLag
    info.terminationCriterionStr = ...
      ['Best-so-far has not been updated for the last gLag generations ', ...
       sprintf('(gLag = %d)', gLag)];
    info.terminationCriterion = 5;
  elseif nFitnessEvaluations + nConstraintEvaluations > config.maxEvals
    info.terminationCriterionStr = ...
      'Maximum number of function evaluations exceeded';
    info.terminationCriterion = 6;
  end
  info.numGenerations = g;
  info.nEvals = nFitnessEvaluations + nConstraintEvaluations;

end

