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

function runConstraintActiveMaEsOnCocoProblem(...
                             problem, lower_bounds, upper_bounds, budget)
  config = struct;
  config.dimension = length(lower_bounds);

  try
    config.xInit = cocoProblemGetInitialSolution(problem);
    config.xInit = config.xInit(:);
    config.maxEvals = budget;
    [x, info] = ConstraintActiveMaEs(...
                                     @(x) fWrapper(problem, x, budget), ...
                                     @(x) gWrapper(problem, x, budget), ...
                                     config);
    %disp('x'), disp(x);
    %disp('info'), disp(info);
  catch
    id = lasterror.identifier;
    if (strcmp(id, 'MAES:budget_exhausted'))
      disp(lasterror.message);
    else
      rethrow(lasterror);
    endif
  end_try_catch
end

function y = fWrapper(problem, x, budget)
  y = cocoEvaluateFunction(problem, x(:));
  if (cocoProblemGetEvaluations(problem) + ...
      cocoProblemGetEvaluationsConstraints(problem) > budget)
    error('MAES:budget_exhausted', 'MA-ES budget exhausted');
  end
end

function y = gWrapper(problem, x, budget)
  y = [];
  if cocoProblemGetNumberOfConstraints(problem) > 0
    y = cocoEvaluateConstraint(problem, x(:));
    if (cocoProblemGetEvaluations(problem) + ...
        cocoProblemGetEvaluationsConstraints(problem) > budget)
      error('MAES:budget_exhausted', 'MA-ES budget exhausted');
    end
  end
  y = y(:);
end

