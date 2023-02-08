% Gaussian computes the values of the Gaussian RBF centered
% around xc with scaling factor scale in the set of points x.
%       Gaussian(x) = psi(||x-xc||/scale)
%            psi(y) = e^(-y^2)
%
% Input:
% - x: (m x d)-array containing the cartesian coordinates of points to
%   evaluate the RBF at
% - xc: (1 x d)-array containing the cartesian coordinates of the center
% - scale (double): scaling factor for the given RBF
%
% Result:
% - (m x 1)- array containing the values of the RBF in x

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function FuncVals = Gaussian(x, xc, scale)

    [npoints, ndim] = size(x);

    % x-xc
    aux = x - repmat(xc, npoints, 1);
    FuncVals = zeros(npoints,1);
    
    %||x-xc||^2 stored in FuncVals
    aux = aux.*aux;
    for i = 1:ndim
        FuncVals = FuncVals + aux(:,i);
    end
    % vectorized version is actually slower than the version above!
    %    FuncVals = aux*ones(ndim,1);
    
    aux2 = 1/(scale*scale);
    FuncVals = exp(-FuncVals*aux2);
end