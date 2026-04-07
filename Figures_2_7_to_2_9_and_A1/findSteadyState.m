function X_star = findSteadyState(params, reduced, X0)
% This function searches for steady state locations in the Thornley model
% by finding points where the right hand side is zero. The user can specify
% whether to naively solve the case where all the RHS equations are zero,
% or the "reduced" case where the first two equations are removed from the
% equations for the following terms. The user provides an initial guess X0,
% which is trialled as an initialisation before further initialisations are
% also trialled in the event a steady state was not found

% Specify the tolerance at which a steady state is considered found
tol = 1e-8;

% Specify the search range to use for random initialisations
X_min = [1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6]';
X_max = [10, 5, 2, 2, 2, 2]';
% Specify the number of initial conditions to try
max_iters = 1000;

% Prepare options to ensure a thorough search
options = optimoptions('fsolve');
options.MaxIterations = 20000;
options.MaxFunctionEvaluations = 20000;
options.OptimalityTolerance = 1e-9;
options.Display = 'off';

% Create a function that is the 2-norm of the right hand side vector, such
% that the minimum occurs where all right hand side elements are zero
if reduced
    min_fun = @(X) reducedRHS(X,params);
else
    min_fun = @(X) thornleyRHS(0,X,params);
end

% Try the search using the default initialisation
X_star = fsolve(min_fun, X0, options);
best_found = norm( min_fun(X_star) );
best_X = X_star;

if best_found <= tol
    searching = false;
else
    searching = true;
end

% Search if a suitable solution has not yet been found
iters = 0;
while searching
   
    % Increment iterations count
    iters = iters + 1;
    % Try the search using a random initial point
    X0 = X_min + rand(6,1) .* (X_max - X_min);
    X_star = fsolve(min_fun, X0, options);
    % Check how close this point gets to zero
    if norm( min_fun(X_star) ) < best_found
        best_found = norm( min_fun(X_star) );
        best_X = X_star;
    end
    
    % Terminate if suitably good point, or too many iterations
    if best_found <= tol || iters >= max_iters
        searching = false;
    end

end

% Return the best point
X_star = best_X;

% Warn user if this point is not sufficiently good
if best_found > tol
    warning(['A sufficiently good steady state was not found, norm of the right hand side vector was ',num2str(best_found)]);
end
    
end