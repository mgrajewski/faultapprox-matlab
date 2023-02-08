% This function tests if a boundary component is closed after adding a
% point to prevent endless loops in prologateFaultLine. We test if any
% point on a remote part of the fault line is in fact close to the
% last point. To do so, we exploit that the points on the fault line
% are already ordered. We take the first or last half of the points on
% the fault line depending on if we prolongate the start or the end
% of the fault line. We compute their distance to the last point. In
% prolongation, this is the point added at last. If this distance is
% small, we conclude that the fault line must be closed.
% This does not make sense for boundary components consisting of
% 6 points only or even less.
%
% Input:
% Points: (npoints x 2)-array of points defining the fault line component
% IdxPoints: index array, Points(IdxPoints) is sorted.
% distVec: array containing the distance between consecutive points in
%   Points
% npoints: number of points in Points
% mode: Mode of operation in prolongation (0: prolongate the current
%   starting point, 1: extend it beyond the current end)
%
% Output:
% closedComp: true, if the fault line component is closed, false, if not

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [closedComp] = testForClosedComp(Points, IdxPoints, distVec, npoints, mode)

    closedComp = false;
    ndim = size(Points, 2);

    if (npoints > 6)
        auxPoints = min(npoints-1, floor(0.5*npoints));

        if (mode == 0)
            auxMax = npoints - 1;
            auxMin = max(1, auxMax - auxPoints);
            distToTestPoints = Points(IdxPoints(auxMin:auxMax) ,:) - Points(end,:);
            distVecTest = distVec(IdxPoints(auxMin:auxMax));

        elseif (mode == 1)
            distToTestPoints = Points(IdxPoints(1:auxPoints) ,:) - Points(end,:);
            distVecTest = distVec(IdxPoints(1:auxPoints));
        end

        % compute the distance of the new points to the test set
        for idim = 1:ndim
            distToTestPoints(:,idim) = distToTestPoints(:,idim).*distToTestPoints(:,idim);
        end

        for idim = 2:ndim
            distToTestPoints(:,1) = distToTestPoints(:,1) + distToTestPoints(:,idim);
        end
        distToTestPoints(:,1) = sqrt(distToTestPoints(:,1));

        % The factor 0.7 has the following reason: if two components really
        % overlap, the minimal distance is at most
        % 0.5*maxDistForSurfacePoints, if a point is on the line segment
        % between two consecutive points. However, as the fault
        % line may be curved, the true minimal distance may be a little
        % larger.
        if (any(distToTestPoints(:,1)' < 0.7*distVecTest))
            closedComp = true;
        end
    end

end
