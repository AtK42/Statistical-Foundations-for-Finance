% This code belongs to the paper
%
% M. Hasannasab, J. Hertrich, F. Laus, and G. Steidl. 
% Alternatives to the EM algorithm for ML-estimation of location, scatter
% matrix and degree of freedom of the student-t distribution.
% Numerical Algorithms, 2020.
% DOI: https://doi.org/10.1007/s11075-020-00959-w
%
% If you use this code, please cite the paper.
%
% This function implements newtons method for determining the zeros in the nu-steps.
%
% INPUTS:
%   init_val        - initial value
%   f               - function handle of the function
%   derivative_f    - function handle of the derivative
%   max_steps       - maximum number of steps
%   tol             - stopping criteria. Stop if abs(f(val))<tol.
%
% OUTPUTS:
%   zero            - root of f
%   evals           - number of evaluations of f and its derivative.
%
function [zero,evals]=ex1a_newton(init_val, f, derivative_f, max_steps, tol)

% set some default values
if nargin < 5
    tol=1e-5;
end
if nargin < 4
    max_steps=1000;
end

% initialize variables
rep=0;
zero=init_val; % note that as initial value, the current nu_r is used
steps=0; 
f_val=f(zero); % apply function on p. 92 (as coded in ex1a_nu_step_mmf.m) on the current nu_r as argument
eps=abs(f_val);

% find zero by using Newton-Raphson (as long as maximum iterations were not exceeded
% and the absolute value of the function value is greater or equal to the
% tolerance as defined above)
while steps<max_steps && eps>=tol
    newzero=zero-f_val/derivative_f(zero); % apply Newton-Raphson
    zero=newzero;
    if zero<0
        zero=10^(-2); % if current potential root is negative (note that we're solving for nu which should be positive), use .01 instead
        rep=rep+1; % increase counter
        if rep>10 % if potential zero is negative 'rep' times, stop the algo
            disp('One nu was reprojected on the interval [1e-2,inf[');
            break;
        end
    end
    f_val=f(zero); % calculation of function value of potential root
    eps=abs(f_val); % take absolute value of function value
    steps=steps+1;
    if mod(steps,100)==0
        disp(['Reached step ' num2str(steps) ' in Newton. Eps: ' num2str(eps) ' with argument ' num2str(zero) ' and derivative ' num2str(derivative_f(zero))])
    end
end

% display message if algo did not converge, i.e., the zero could not be found within the set bounds
if eps>=tol
    disp(['Newton did not converge in ' num2str(max_steps) ' steps. Eps: ' num2str(eps) ' with argument ' num2str(zero) ' and derivative ' num2str(derivative_f(zero))])
end

% additionally to the zero, give out how many iterations were needed
evals=2*steps+1;

disp(zero);
end