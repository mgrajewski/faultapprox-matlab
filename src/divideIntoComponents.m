% This function separates a given fault line in its components. To do so,
% it relies on the fact that PointSet is ordered. This functions applies to
% 2D only. The function detects several components, if two subsequent
% points have a mutual distance greater than maxDistForSurfacePoints. Note
% that this function is called AFTER algorithm fill.
%
% Input:
% - PointSet: set of points representing the fault
% - maxDistForSurfacePoints: maximum distance of two subsequent triplets on
%   a fault line
%
% Output:
% - Components: cell arrays of arrays containing representing sets of the
%   components
% - numOfComps: the number of components the fault line consists of
% - NumPointsPerComp: array containing the number of points per component

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [Components, numOfComps, NumPointsPerComp] = ...
    divideIntoComponents(PointSet, maxDistForSurfacePoints)

    numPoints = size(PointSet,1);
    dist = PointSet(1:numPoints-1,:) - PointSet(2:numPoints,:);
    dist(:,1) = dist(:,1).^2 + dist(:,2).^2;
    dist(:,1) = sqrt(dist(:,1));
    CompAux = dist(:,1) > 3*maxDistForSurfacePoints;
    numOfComps = size(CompAux(CompAux > 0),1)+1;

    NumPointsPerComp = zeros(1, numOfComps);
    compEnd = find(CompAux > 0);

    Components = cell(numOfComps, 1);

    compStart = 1;
    for icomp = 1: numOfComps-1
        Components{icomp} = PointSet(compStart: compEnd(icomp),:);
        NumPointsPerComp(icomp) = compEnd(icomp) - compStart + 1;
        compStart = compEnd(icomp) + 1;
    end
    Components{numOfComps} = PointSet(compStart:numPoints,:);
    NumPointsPerComp(numOfComps) = numPoints - compStart + 1;
end