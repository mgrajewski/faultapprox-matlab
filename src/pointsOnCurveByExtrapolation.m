% This function computes numPointsNew points on an extrapolated fault
% line. It is constructed using the points given in PointsOnCurve and
% fitting a polynomial in local coordinates. As the points are not
% known exactly, we do not interpolate exactly, but penalise the second
% derivative in a least-squares approximation. Following Morozov, we
% choose the regularisation such that the maximal residual is approx.
% FaultApproxParams.abstolBisection.
% These points are constructed such they have approximately the same 
% distance as the points on PointsOnCurve, but at most
% maxDistForSurfacePoints.
%
% Input:
% - PointsOnCurve: points on a curve (must be ordered)
% - numPointsNew: number of new points on the curve aka fault line
%   constructed by extrapolation.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - NewPointsOnCurve: cartesian coordinates of the new points
% - avgDist: average distance of the points on the fault line used for
%   extrapolation
% - DataPars: Parameters aka local coordinates of points in PointsOnCurve
% - ExtraPars: Parameters aka local coordinates of the new points on the
%   polynomial curve
% - Q: Orthogonal matrix; the local coordinates are computed by
%   xloc = Q*x + xmean
% - xmean: origin of the local coordinate system; the local coordinates
%   are computed by xloc = Q*x + xmean
% - coeffs: coefficients of the extrapolating polynomial

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [NewPointsOnCurve, avgDist, ExtraPars, DataPars, Q, xmean, coeffs] = ...
    pointsOnCurveByExtrapolation(PointsOnCurve, numPointsNew, FaultApproxParams)
    
    growthFactor = 1.5;
    
    [numPointsOnCurve, ndim] = size(PointsOnCurve);

    % we need at least two points for extrapolation
    if (numPointsOnCurve < 2)
        warning('Too few points to extrapolate')
        NewPointsOnCurve = [];
        avgDist = 0;
        DataPars = [];
        Q = zeros(ndim);
        xmean = PointsOnCurve;
        coeffs = [];
        return
    end
    
    % average distance of the points on the fault line used for
    % extrapolation
    avgDist = norm(PointsOnCurve(numPointsOnCurve,:) - PointsOnCurve(1,:),2)/(numPointsOnCurve-1);
    
    % mean value of coordinates
    xmean = 1/size(PointsOnCurve, 1)*ones(1, size(PointsOnCurve, 1))*PointsOnCurve;

    % shift point set such that its mean is zero                   
    PointsShifted = PointsOnCurve - xmean;
    [~,~,Q] = svd(PointsShifted);
    
    % rotate such that the points are aligned to the x-axis
    PointsRot = PointsShifted*Q;

    % We limit the stepsize for interpolation based upon the curvature. It
    % does not make sense to extrapolate very far with a lot of numerical
    % amount and to adaptively refine in between the current and the new
    % points anyway afterwards. To do so, we compute the step size which
    % would lead to the maximal admissible deviation
    % FaultApproxParams.errMax from a straight line segment. This is a
    % natural upper bound for the step size in extrapolation.

    % The maximal deviation d of a curve with curvature c from an
    % interpolating line segment with length l is
    % d = 1/4 c l^2 + 1/16 c^3 l^4 + h.o.t
    % Rearranging leads to l = 2/c^2(sqrt{1+4cd} - 1), which is numerically
    % unstable if cd is small. We set (cl)^2 = v and search for the bigger
    % root of the quadratic equation
    % v^2 + 4 v - 16cd = 0, aka v^2 + 2pv - q = 0. According to
    % Vieta, q = vmin*vmax with the two roots vmin and vmax. Therefore,
    % vmin = - (2 + sqrt{4+ 16cd}) = -2 (1 + sqrt{1 + 4cd}), which is stable
    % to evaluate and vmax = -16cd/vmin, and ultimately
    % lmax = sqrt{vmax}/c. Inserting vmax yields
    % lmax = sqrt{-16cd/vmin}/c = 4 sqrt{d/(-c vmin)}
        
    % estimating curvature is possible only for more than two points on the
    % fault line and two dimensions
    if (numPointsOnCurve > 2 && ndim == 2)

        % estimate curvature (with safety factor). As PointsOnCurve are mean
        % values of the point pairs, their maximal deviation from the true
        % fault line is FaultApproxParams.abstolBisection.
        curv = computeCurvature(PointsShifted, numPointsOnCurve, FaultApproxParams.abstolBisection) + ...
               computeCurvature(PointsShifted, numPointsOnCurve-1, FaultApproxParams.abstolBisection);

        if (curv > 1e-10)
            lmax = 2*(1 + sqrt(1 + 4*curv*FaultApproxParams.errMax));
            lmax = 4*sqrt(FaultApproxParams.errMax/(curv*lmax));
        else
            lmax = 1e10;
        end
        
    % limiting step size based upon curvature does not work, if only two
    % points on the fault line are known    
    else
        lmax = 1e10;
    end
           
    stepSize = min([FaultApproxParams.maxDistForSurfacePoints, growthFactor*avgDist, lmax]);
    
    if (all(PointsRot(1:numPointsOnCurve-1,1) < PointsRot(2:numPointsOnCurve,1)))    
        ExtraPars = PointsRot(numPointsOnCurve,1) + stepSize * [1:numPointsNew]';
    elseif (all(PointsRot(1:numPointsOnCurve-1,1) > PointsRot(2:numPointsOnCurve,1)))
        ExtraPars = PointsRot(numPointsOnCurve,1) - stepSize * [1:numPointsNew]';

    % It is not possible to extrapolate the fault line by any function
    % respecting the order of the points. This is usually the case if the
    % points are on a curve with very strong curvature. Therefore, we just
    % take the two last points and perform linear extrapolation.
    else
        numPointsOnCurve = 2;
        if (PointsRot(numPointsOnCurve-1,1) < PointsRot(numPointsOnCurve,1))
            ExtraPars = PointsRot(numPointsOnCurve,1) + stepSize*[1:numPointsNew]';
        else
            ExtraPars = PointsRot(numPointsOnCurve,1) - stepSize*[1:numPointsNew]';
        end
    end
    
    A = zeros(numPointsOnCurve);
    B = zeros(numPointsOnCurve);
    for i = 1: numPointsOnCurve
        A(:,i) = PointsRot(1:numPointsOnCurve,1).^(i-1);
        B(:,i) = (i-1)*(i-2)*PointsRot(1:numPointsOnCurve,1).^(max(0,i-3));
    end
    
    % Regularization by penalizing the second derivative is pointless if
    % only two points exist so far: Then, we perform linear extrapolation,
    % and its second derivative is zero anyway.
    if (numPointsOnCurve > 2)
    
        ATA = A'*A;
        BTB = B'*B;
        rhs = A'*PointsRot(:,2);

            % estimation of regularisation parameter due to Morozov
        expmin = -16;
        expmax = 2;

        eps = FaultApproxParams.abstolBisection;
        nit = 0;
        while (expmax - expmin > 1/3)

            nit = nit+1;

            expnew = 0.5*(expmin + expmax);
            munew = 10^expnew;

            Awork = ATA;
            Awork = Awork + munew*BTB;

            coeffs = linsolve(Awork, rhs);

            res = norm(A*coeffs - PointsRot(:,2));

            if (res > eps)
                expmax = expnew;
            else
                expmin = expnew;
            end

            % small enough
            if (munew < 1e-14)
                break
            end
        end
    else
        coeffs = linsolve(A,PointsRot(1:numPointsOnCurve,2));
    end
    
    NewPointsOnCurve = evalPol(ExtraPars, coeffs, Q, xmean);
    DataPars = PointsRot(:,1);
    
    trueDist = vecnorm(NewPointsOnCurve- PointsOnCurve(numPointsOnCurve,:),2,2);
    
    % it may happen that the true distance of the extrapolated points to
    % the existing ones is too large, as the arc length on the curve is
    % much greater than the distance of the parameters. In this case, we
    % adjust the parameters accordingly. For the sake of simplicity, we
    % consider the true distance of the last existing to the first new
    % point only.
    if (trueDist(1) > 1.1*FaultApproxParams.maxDistForSurfacePoints)

        tmax = ExtraPars(1);        
        tmin = DataPars(numPointsOnCurve);

        iiter = 1;

        % the maximal number of iterations should never be reached in
        % practical examples. If so, the adjusted values are at least
        % better than the original ones.
        while iiter < 20
            tnew = (tmin+tmax)/2;
            xnew = evalPol(tnew, coeffs, Q, xmean);
            
            dist = norm(xnew - PointsOnCurve(numPointsOnCurve,:));
            res = dist - FaultApproxParams.maxDistForSurfacePoints;
            
            % if the distance to the existing points is
            % maxDistForSurfacePoints up to 10%
            if (abs(res) < 0.1*FaultApproxParams.maxDistForSurfacePoints)
                break
            else
                if dist > FaultApproxParams.maxDistForSurfacePoints
                    tmax = tnew;
                else
                    tmin = tnew;
                end
            end
            iiter = iiter + 1;
        end
       
        ExtraPars = DataPars(numPointsOnCurve) + [1:numPointsNew]'*(tnew-DataPars(numPointsOnCurve));
        NewPointsOnCurve = evalPol(ExtraPars, coeffs, Q, xmean);
    end
end