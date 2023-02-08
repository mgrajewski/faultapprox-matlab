% We assume that a planar curve is locally represented by a set of points,
% PointSet. This function estimates the curvature in the point in PointSet
% with index iidx. The point set is assumed to be ordered, however, the
% points do not need to be equidistant.
% After shifting and rotation, we assume that the curve can be represented
% as a graph of an unknown function. It is straightforward to interpolate
% this graph using Gaussian Radial Basis Functions and to estimate the true
% curvature by the one of this interpolating function.
% In the context of fault detection, however, the points in PointSet are
% located on the curve only up to a tolerance epsTol. Therefore, we do not
% interpolate, but penalise the second derivative of the interpolating RBF
% function subject to a maximal residual < epsTol. Note that the maximal
% residual coincides with the maximal deviation in the value at an
% interpolation point. As the points are known up to epsTol, it is 
% pointless to interpolate more exactly. For computing this, we employ
% Tikhonov regularization with parameter estimation following Morozov.
%
% If the curve cannot be considered a graph even after rotation, we draw as
% a fallback a circle through the points with indices iidx-1, iidx and 
% iidx+1 and use its radius for estimating the curvature.
% It may happen that in case of extremely unevely distributed points, the
% unregularized matrix A^TA is ill conditioned. Unevely distributed points 
% may occur, if the fault line has a kink, adaptive refinement/coarsening
% concentrates the points at the kink and eliminates points far from it.
% However, the regularisation with the second derivative reduces the 
% condition number of the interpolation matrix to reasonable size, such 
% that we can ignore the warnings issued by the linear algebra library.
%
% For efficiently inserting new points during adaptive refinement, it
% makes sense to construct starting pairs based on the RBF approximation
% at hand. Therefore, we provide the (approximate) left and right midpoints
% from the point with index on the approximating curve, if such points
% exist. These points are given in global coordinates.
%
% Input:
% - PointSet: Set of points (size npoints x ndim)
% - iidx: index in PointSet of the point to estimate the curvature at
% - epsTol: tolerance for the residual in regularisation
%
% Output:
% - curvature: estimated curvature in point iidx
% - NewPoints: coordinates of the left and right midpoints of iidx on
%   the RBF curve.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [curvature, NewPoints] = computeCurvature(PointSet, iidx, epsTol)
    maxIt = 30;
    NewPoints = 0;

    [numPoints, ndim] = size(PointSet);

    if (numPoints < 3)
        message = ['Not enough points (at least three needed) ', ...
                   'given for estimating the curvature'];
        error(message)
    end
    
    if (ndim ~= 2)
        error('computeCurvature works in 2d only.')
    end
    
    if (iidx < 1 || iidx > numPoints)
        message = ['The index of the point to estimate the curvature ' ...
                   'at is either  below 0 or larger than npoints.'];
        error(message)
    end
    
    % first step: local transformation
    
    % shift such that the point to estimate the curvature at is 0
    oldCenter = PointSet(iidx,:);
    PointSet = PointSet - oldCenter;

    % Compute local coordinate system based on the principal axes
     
    % Orthogonal transformation: we consider the rightmost and leftmost
    % point. Let alpha be the angle between the line connecting these
    % points and the x-axis. We turn all points in PointSet by -alpha
    % degrees applying Q with c = cos(alpha) and s = sin(alpha).
    % Of course, this is a heuristic approach which can fail.
    cs = (PointSet(end,:)- PointSet(1,:))/ ...
          norm(PointSet(end,:) - PointSet(1,:));
    Q = [cs(1), -cs(2); cs(2), cs(1)];
    
    % Transform the point set (note that rotations preserve curvature)
    % to be aligned to the x-axis as good as possible. Therefore, the
    % problem to find the curvature of a curve has been reduced to compute
    % the curvature of a graph.
    PointSet = PointSet*Q;
    
    % We assume that, after rotation, the part of the line can be
    % considered as a graph of a function. As the points have been
    % presorted, the x-values must be descending or ascending after
    % transformation.
    if ~(all(PointSet(1:numPoints-1,1) < PointSet(2:numPoints,1)) || ...
         all(PointSet(1:numPoints-1,1) > PointSet(2:numPoints,1)))
        % Obviously, curvature is very strong in this case. So, we draw a
        % circle through three consecutive points and use its radius for
        % estimating curvature. This situation occurs however for
        % consecutive points only if the resolution is far too coarse for
        % the fault line or if the sorting is wrong. In any case, we need
        % more points in this region. Therefore, we just return the value
        % of the curvature, which should be rather large even if wrong to
        % enforce further refinement in this region.
        if (iidx < numPoints && iidx > 1)
            %a = PointSet(iidx-1, :);
            %b = PointSet(iidx+1, :);

            [curvature,~] = curvCenterFromThreePoints(PointSet(iidx-1: iidx+1, :));
            NewPoints = 0.5*[PointSet(iidx-1, :); PointSet(iidx+1, :)];

        elseif (iidx == 1)
            % center PointSet around second point
            auxCenter = PointSet(2,:);
            PointSet = PointSet - auxCenter;

            [curvature,~] = curvCenterFromThreePoints(PointSet(1: 3, :));
            NewPoints = 0.5*PointSet(1, :) + auxCenter;

        elseif (iidx == numPoints)

            % center PointSet around second to last point
            auxCenter = PointSet(iidx-1,:);
            PointSet = PointSet - auxCenter;

            [curvature,~] = curvCenterFromThreePoints(PointSet(numPoints-2:numPoints, :));
            NewPoints = 0.5*PointSet(iidx-2, :) + auxCenter;
        end
        NewPoints = NewPoints*Q' + oldCenter;
        
        return;        
    end
    
    % build interpolation matrix
    % It depends on the x-values of the rotated points only, the y-values
    % contribute to the right hand side
    A = zeros(numPoints, numPoints);
    ScaleVec = abs((PointSet(1,1) - PointSet(end,1)))/ ...
               numPoints*ones(numPoints, 1);
    
    for i = 1:numPoints
        A(:,i)= Gaussian(PointSet(:,1), PointSet(i,1), ScaleVec(i));
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
    B = zeros(numPoints, numPoints);
    phi = zeros(1,numPoints);
    for i = 1:numPoints
        B(:,i) = Gaussian_second_der(PointSet(:,1), PointSet(i,1), ...
                                     ScaleVec(i));
    end
    B = B'*B;
    
    % if the points are close together, B is very large, albeit well
    % conditioned. We normalise B in order to balance the size of A and B.
    B = 10/norm(B)*B;
    
    % fit RBF curve
    ATA = A'*A;  
    rhs = A'*PointSet(:,2);

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
        res = norm(A*coeffs - PointSet(:,2), inf);

        if (res > epsTol)
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
    
    % x-positions of the new points: arithmetic mean between current
    % points and its left and right neighbour, if it exists. Note that the
    % coordinates of the point with index iidx are shifted to 0.
    if (iidx > 1 && iidx < numPoints)

        auxPos = 0.5*[PointSet(iidx-1,1); PointSet(iidx+1,1)];
    
    elseif (iidx == 1)
        auxPos = 0.5*PointSet(iidx+1,1);
    elseif (iidx == numPoints)
        auxPos = 0.5*PointSet(iidx-1,1);        
    end

    phival = zeros(size(auxPos,1), numPoints);
    for i = 1: numPoints
        phival(:, i) = Gaussian(auxPos, PointSet(i,1), ScaleVec(i));
    end
    NewPoints = phival*coeffs;
    NewPoints = [auxPos, NewPoints]*Q' + oldCenter;

    % compute first derivative
    for i = 1: numPoints
        phi(i) = Gaussian_first_der(PointSet(iidx,1), PointSet(i,1), ...
                                    ScaleVec(i));
    end
    phidot = phi*coeffs;
    
    % compute second derivative
    for i = 1: numPoints
        phi(i) = Gaussian_second_der(PointSet(iidx,1), PointSet(i,1), ...
                                     ScaleVec(i));
    end
    phiddot = phi*coeffs;
    
    curvature = abs(phiddot/(1 + phidot^2)^1.5);
end