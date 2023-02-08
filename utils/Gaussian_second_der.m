% Gaussian_second_der computes the values of the second derivative of the
% Gaussian RBF centered around xc with scaling factor scale in the set of
% points x.
%       Gaussian(x) = psi(||x-xc||/scale)
%            psi(y) = e^(-y^2)
%
%       Gaussian_second_der(x) =
%       psi(||x-xc||/scale)*( 4/(scale^4)*(x-xc)*(x-xc)' - 2/(scale^2)*I)
%
% Input:
% - x: (m x d)-array containing the cartesian coordinates of points to
%   evaluate the RBF at
% - xc: (1 x d)-array containing the cartesian coordinates of the center
% - scale (double): scaling factor for the given RBF
%
% Result:
% - (m x d x d)- array containing the values of the second derivative of
%   the RBF in x

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function FuncVals = Gaussian_second_der(x, xc, scale)

    [npoints, ndim] = size(x);

    % x-xc
    aux = x - repmat(xc, npoints, 1);
    
    %||x-xc||^2 stored in auxNorm2
    auxNorm2 = zeros(npoints,1);
    auxNorm = aux.*aux;
    for i = 1:ndim
        auxNorm2 = auxNorm2 + auxNorm(:,i);
    end
    
    FuncVals = zeros(npoints, ndim, ndim);
    aux2 = 1/(scale*scale);
    
    for ipoint = 1: npoints
        FuncVals(ipoint,:,:) = (4*aux2*aux2*(aux(ipoint,:)'*aux(ipoint,:)) - 2*aux2*eye(ndim))*exp(-auxNorm2(ipoint,:)*aux2);
    end
end