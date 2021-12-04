% FUNCTION NAME:
%   HorCyl_AreaFromLevel
%
% DESCRIPTION:
%   Compute the cross-sectional area of a horizontal, cylindrical tank,
%		corresponding to a given level (height)
%
% INPUT:
%   h   - (double) Height of filled segment (level)     [m]
%   D   - (double) Diameter of circular cross-section	[m]
%
% OUTPUT:
%   A   - (double) Calculated cross-sectional area corresponding to level
%
function A = HorCyl_AreaFromLevel( h, D )
  
  % Calculate area
  f = ( asin(2*h/D - 1) + 2*(2*h/D - 1)*sqrt(h/D - h*h/D/D) ) / pi + 1/2;
  A = f * pi * D*D/4;

end
