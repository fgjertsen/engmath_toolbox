% FUNCTION NAME:
%   nleq_broyden
%
% DESCRIPTION:
%   Find the solution for a nonlinear (set of) equation(s)
%		on the form Ax = b, by iteration
%
% INPUT:
%   xv0	- (double, vector) 				Initial guess for solution
%   f   - (fhandle(double, vector))		Function to solve, f(x, par) = 0
%	par	- (double, struct)				System specific parameters
%	J 	- (fhandle(double, matrix)) 	Function to calculate Jacobian, df/dx.	(optional)
%	tol	- (double, scalar)				Tolerance for iteration convergence 	(optional)
%
% OUTPUT:
%   xv	- (double, vector) 	Solution to equation
% 	nIt - (integer, scalar)	Number of iterations needed to converge
%
function [xv, nIt] = nleq_broyden(xv0, f, par, J, tol)

	% Initiate function values
	xv = x0;
	n  = length(xv);
	fr = feval(f, xv, par);
	
	% Initiate Broyden approximation to Jacobian
	if not( isempty ( J ) )
		Br = inv( feval( J, xv, par) );
	else
		Br = eye(n);
	end
	if isempty(tol)
		tol = 1.e-10;
	end
	M = zeros(n);
	
	% Tuning parameter of Broydens method
	tau = 1;
	
	% Start iteration loop
	bConverged = 0;
	nIt = 0;
	while ~conv
		
		% Check for convergence
		if norm(fr) < tol
			bConverged = 1;
			break;
		end
		
		% Do "Newton Iteration" step
		pr = -Br * fr;
		xv = xv + tau*pr;
		
		% Update Jacobian approximation, using Broydens formula
		dfr = feval(f, xv, par) - fr;
		fr  = fr + dfr;
		oyp = Br*y-pr;
		pB	= pr' * Br;
		for i=1:n
			for j=1:n
				M(i,j) = oyp(i) * pB(j);
			end
		end
		Br = Br - M./(pr'*Br*dfr);
		
		% Increment loop counter
		nIt = nIt + 1;
		
	end

end