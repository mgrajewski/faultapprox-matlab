% This function aims at finding the intersection of an extrapolated fault
% line with the domain boundary using bisection. We exploit that the domain
% is either a rectangle or a cuboid.
% We assume that the fault line has been extrapolated with a polynomial in
% local coordinates. We approximate the parameter value t of the
% intersection of that polynomial with the domain boundary by bisection.
% For doing so, we need a lower and upper bound tmin and tmax for that
% parameter t along with the corresponding points xmin and xmax.
%
% Input:
% - tmin, tmax: lower and upper bound for the parameter of the intersection
% - xmin, xmax: corresponding points in global coordinates
% - coeffs: coefficients of the extrapolating polynomial
% - Q, xmean: orthogonal matrix and origin describing the local coordinate
%   system
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - xnew: approximate intersection of fault line and domain boundary
% - iedge: index of the domain boundary part, at which the fault line
%   leaves the domain

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [xnew, iedge] = extraNewPointOnDomainBdry(tmin, tmax, xmin, xmax, ...
                                                   coeffs, Q, xmean, ...
                                                   ProblemDescr, ...
                                                   FaultApproxParams)

    % indices of the domain edges (in 2D):
    %        y /\
    %          |       3
    %          |_______________ 
    %          |               |
    %          |               |
    %          |               |
    %        4 |               | 2
    %          |               |
    %          |               |
    %          |_______________|__________ x
    %                  1
    %
    % indices of the domain facets (in 3D)
    %
    %                  6
    %                  |     3                  z
    %             _____|____/_____              |   y
    %            /|    |   /     /|             |  /
    %           / |    V  /     / |             | /
    %          /  |     |/_    /  |             |/
    %         /   |           /   |             +---------x
    %        /_______________/    | 
    %        |    |          |  <--------2
    %  4-----|->  |__________|____|
    %        |   / _         |   /
    %        |  /  /|  /\    |  /
    %        | /  /    |     | /
    %        |/__/_____|_____|/
    %           /      |
    %          1       5 
    iauxMin = [4 1 5];
    iauxMax = [2 3 6];
    
    % points closer as epstol are regarded as identical
    epstol = FaultApproxParams.eps;

    ndim = size(xmin,2);

    % euclidean distance of these points
    dist = norm(xmin - xmax);
    iiter = 1;
    while iiter < 10

        % the 0.5 is just a safety factor (iterating one time more does
        % not increase the number of function evaluations)
        if (dist < 0.5*FaultApproxParams.abstolBisection)
            break
        else
            tnew = 0.5*(tmin + tmax);                    
            xnew = evalPol(tnew, coeffs, Q, xmean);
        
            if (any(xnew < ProblemDescr.Xmin) || ...
                any(xnew > ProblemDescr.Xmax))
                tmax = tnew;
                xmax = xnew;
            else
                tmin = tnew;
                xmin = xnew;
            end
            
            % The distance will approximately halven per iteration step,
            % but only approximately (the extrapolation is not necessarily
            % a straight line).
            dist = norm(xmax - xmin);
        end
        iiter = iiter + 1;
    end

    % index of the edge or facet at which the fault line leaves
    % the domain
    iedge = 0;
    for i = 1:ndim
        if (xmax(1,i) < ProblemDescr. Xmin(i))
            iedge = iauxMin(i);
        elseif (xmax(1,i) > ProblemDescr. Xmax(i))
            iedge = iauxMax(i);
        end
    end
    
    % The approximation of the intersection is not necessarily
    % inside the domain, so we enforce this here.
    xnew = max(min(xnew, ProblemDescr.Xmax - epstol), ...
               ProblemDescr.Xmin + epstol);

    % Moreover, ensure that the new point is up to epstol on the
    % domain boundary (it is somewhere inside near a boundary
    % up to abstolBisection, but we can do better).
    switch iedge
        
        % x2 is minimal (bottom boundary/front surface of cube)
        case 1
        xnew(2) = ProblemDescr.Xmin(2) + epstol;

        % x1 is maximal (right boundary/right surface of cube)
        case 2
        xnew(1) = ProblemDescr.Xmax(1) - epstol;

        % x2 is maximal (top boundary/back surface of cube)
        case 3
        xnew(2) = ProblemDescr.Xmax(2) - epstol;

        % x1 is minimal (left boundary/left surface of cube)
        case 4
        xnew(1) = ProblemDescr.Xmin(1) + epstol;

        % x3 is minimal (bottom surface of cube)
        case 5
        xnew(3) = ProblemDescr.Xmin(2) + epstol;

        % x3 is maximal (top surface of cube)
        case 6
        xnew(3) = ProblemDescr.Xmax(2) - epstol;
    end
end