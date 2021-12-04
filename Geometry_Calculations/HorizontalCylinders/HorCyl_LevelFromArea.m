% FUNCTION NAME:
%   HorCyl_LevelFromArea
%
% DESCRIPTION:
%   Compute the (liquid) level of a horizontal, cylindrical tank, from
%		a given cross-sectional area
%
% INPUT:
%   A   - (double) Height (level)                       [m]
%   D   - (double) Diameter of circular cross-section 	[m]
%
% OUTPUT:
%   h   - (double) Calculated cross-sectional area corresponding to level
%
function L = HorCyl_LevelFromArea( A, D )

	% Set independent variable
	R = D/2;
	h = -R:0.001:R;	
	
	a = A / R / R;
	
	% Establish Newton iterations
	x = 2*a/pi;
	nMaxIter = 50;
	for i = 1:nMaxIter
		fx = acos(1-x) + sqrt(x*(2-x))*(x-1) - a;
		if (abs(fx) < 1.e-30)
			break;
		end
		
		dfx = 2 * sqrt(2*x - x*x);
		
		x = x - fx/dfx;	
	end
	
	inds = find(h < h(1) + x*R);
	
	L = (h(inds(end)) - h(1));

end
