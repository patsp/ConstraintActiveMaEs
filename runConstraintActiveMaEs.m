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

function runConstraintActiveMaEs
  config.dimension = 3;
  config.xInit = 10 * ones(config.dimension, 1);
  config.maxEvals = 100000 * config.dimension;
  [x, info] = ...
  ConstraintActiveMaEs(@(x) objectiveFun(x), ...
                       @(x) constraintFun(x), ...
                       config);
  disp('info'); disp(info);
end

function y = objectiveFun(x)
  y = x' * x;
end

function y = constraintFun(x)
  % x(1) >= 1 => -x(1) + 1 <= 0
  y = [-x(1) + 1];
end

