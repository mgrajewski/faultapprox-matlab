% This function evaluates a polynomial given by the its coefficient coeffs
% in t. The result x is computed in local coordinates defined by Q
% (rotation) and xmean (translation). We however return the result (t,x) in
% global coordinates.
%
% Input:
% - t: parameter value to evaluate the polynomial at
% - coeffs: coefficient vector defining the polynomial
% - Q: orthogonal matrix for the local coordinates
% - xmean: origin of the local coordinates
%
% Output:
% - x: (t, p(t)) in global coordinates

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function x = evalPol(t, coeffs, Q, xmean)
       
    ndim = size(xmean,2);
    npoints = size(t, 1);
    
    deg = size(coeffs,1) - 1;
    x = (t.^(0:deg))*coeffs;

    if (ndim == 2)
        x = [t x]*Q' + xmean;
    else
        x = [t x zeros(npoints,1)]*Q' + xmean;
    end
end