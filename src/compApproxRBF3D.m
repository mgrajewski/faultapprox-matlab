% This function computes an RBF approximation to a function given by its
% values PointValues in PointsInPlane. The values are up to eps_res
% accurate such that we do not just interpolate, but employ Tikhonov
% regularization with parameter estimation due to Morozov.
%
% Input:
% - PointsInPlane: coordinates of the evaluation points
% - PointValues: approximate function values in PointsInPlane
% - eps_res: accuracy of the function values
%
% Output:
% - coeffs: coefficient vector of the RBF approximation
% - scale: shape parameter of the RBF approximation

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [coeffs, scale] = compApproxRBF3D(PointsInPlane, PointValues, eps_res)
    
    maxIt = 30;

    numPoints = size(PointsInPlane, 1);
    
    % fit Gaussian RBFs
    A = zeros(numPoints, numPoints);
    scale = max(PointsInPlane) - min(PointsInPlane);
    scale = 2*sqrt(scale(1)*scale(2)/numPoints);
    ScaleVec = scale*ones(numPoints, 1);
    
    for i = 1:numPoints
        A(:,i)= Gaussian(PointsInPlane, PointsInPlane(i,:), ScaleVec(i));
    end
    
    % penalty matrix
    % we want to penalise the second derivative of the almost interpolating
    % RBF function f, aka ||f''||^2. We approximate this by
    % \sum (f''(x_i))^2, where x_i denote the interpolation points. As
    % f''(x_i) = \sum coeff_j phi_j''(x_i), we can express the vector
    % (f''(x_1), ..., f''(x_n)) by
    %
    % /f''(x_1)\   /phi_1''(x_1) ... phi_n''(x_1) \   /coeff_1\
    % |    .   |   |      .                .      |   |   .   |
    % |    .   | = |      .                .      | * |   .   |
    % |    .   |   |      .                .      |   |   .   |
    % \f''(x_n)/   \phi_1''(x_n) ... phi_n''(x_n) /   \coeff_n/
    %                    = B                            = x
    %
    % Therefore, we have ||f''||^2 \approx <Bx, Bx> = x^T (B^TB) x, such
    % that penalty matrix is B^T B.
    
    % For some unknown reason, penalizing with the 2nd derivative does not
    % work properly here. Therefore, we just take the unit matrix. We
    % suspect a bug.
    B = eye(numPoints);


    % fit RBF curve
    ATA = A'*A;  
    rhs = A'*PointValues;

    % estimation of regularisation parameter due to Morozov
    expmin = -16;
    expmax = 2;
    
    nit = 0;
    while (expmax - expmin > 1/3)

        nit = nit+1;
        
        expnew = 0.5*(expmin + expmax);
        munew = 10^expnew;

        Awork = ATA;
        Awork = Awork + munew*B;

        coeffs = linsolve(Awork, rhs);

        % res is the maximal absolute deviation from the interpolation value
        res = norm(A*coeffs - PointValues, inf);

        if (res > 2*eps_res)
            expmax = expnew;
        else
            expmin = expnew;
        end

        % small enough
        if (munew < 1e-14)
            break
        end
        
        if (nit > maxIt)
            message = ['maximal number of iterations reached for ' ...
                       'searching the regularization parameter'];
            warning(message)
        end
    end        
end