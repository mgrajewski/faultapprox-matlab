% It may happen that between two consecutive triplets on a fault line,
% there is a substantial gap. However, the maximal distance between
% consecutive triplets on a fault line is prescribed by the user
% setting the parameter FaultApproxParams.maxDistForSurfacePoints.
% Therefore, we need to fill such gaps. This is what this function does.
% We sort the given triplets (aka point pairs) on the fault line and 
% compute the distance between consecutive points. If we detect a gap,
% we create new equidistant point pairs. The starting points for
% bisection are computed by the estimated local curvature of the fault
% line. For creating them, we divide the straight line from the current
% existing point on the fault line to its successor equidistantly and
% deviate from it in normal direction symmetrically in both directions.
% The distance to the line is driven by the estimated curvature and
% safeguarded such that for a straight line, no further function
% evaluations apart of the ones for classification of the starting points
% should be necessary. This is for efficiency reasons.
% In the case of strong curvature, it may happen that by coincidence, a
% newly added point pair is very close to an existing one. Before adding a
% new point, we consider the minimal distance to the existing ones and add
% the point only if this distance is greater than
% maxDistForSurfacePoints*minDistFactor.
% We do not resort the points on the line after filling gaps but just
% append them to the arrays of existing points.
%
% Input:
% - PointsIclass, PointsJclass: set of triplets describing the fault line 
% - iclass, jclass: class indices
% - ClassVals: Array containing the class values. These are not
%   necessarily the class indices. Imagine that f(\Omega) = {1,2,5}. Then,
%   ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
%   Size: nclasses
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - NormalVecOfPlane: Outer normal vector of the plane filling takes place
%   in. We need this information, as we reuse this function in expanding
%   fault surfaces in 3D: We treat the intersection of a fault surface with
%   a facet of the domain cuboid as a 2D-fault line and apply this function
%   to that case. 
%
% Output:
% - PointsIclass, PointsJclass: set of triplets describing the fault line
%   (this time without gaps)
% - nTriplets: number of triplets in PointsIclass/PointsJclass
% - tripletsAdded: true, if some points/triplets have been added; false
%   otherwise
% - success: true, if filling of gaps was successful, false otherwise

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointsIclass, PointsJclass, nTriplets, tripletsAdded, success] = ...
    fill2D(PointsIclass, PointsJclass, iclass, jclass, ClassVals, ...
           FaultApproxParams, ProblemDescr, NormalVecOfPlane)

    global ncalls

    safetyFactor = 1.5;
    
    % number of triplets on the fault line
    [nTriplets, ndim] = size(PointsIclass);
    
    maxDistForSurfacePoints = FaultApproxParams.maxDistForSurfacePoints;
    epsSafemax = FaultApproxParams.alpha;
    abstolBisection = FaultApproxParams.abstolBisection;

    AuxVec = zeros(1,3);
    
    % true, if FillGapsInFaultLine added some triplets
    tripletsAdded = false;
    
    [IdxPointsSurfOrdered, sortSuccess] = sortPointsOnFaultLine(PointsIclass, ...
                                          true, ProblemDescr, FaultApproxParams);

    if (~sortSuccess)
        warning(['FillGapsInFaultLine failed, iclass = ' int2str(iclass) ...
                 ', jclass = ' int2str(jclass)])
        success = false;
        return
    end
    
    % resort point sets
    PointsIclass = PointsIclass(IdxPointsSurfOrdered{1},:);
    PointsJclass = PointsJclass(IdxPointsSurfOrdered{1},:);
    
    % here, the points on the fault line are ordered. Fill gaps, if
    % necessary.
    iidxNextSegment = 1;
    IdxSegment = zeros(nTriplets, 1);
    
    StartPointsLeft = zeros(2*nTriplets-1, ndim);
    StartPointsRight = zeros(2*nTriplets-1, ndim);

    % distance to the next point on the fault line
    distNext = vecnorm(PointsIclass(1:nTriplets-1,:) - ...
                       PointsIclass(2:end,:),2,2);
    
    for ipoint = 1: nTriplets - 1
        if(distNext(ipoint) > maxDistForSurfacePoints)

            % number of points to fill the gap
            numSubdivisions = ceil(distNext(ipoint)/maxDistForSurfacePoints);

            IdxSegment(ipoint) = iidxNextSegment;
            iidxNextSegment = iidxNextSegment + numSubdivisions-1;

            parameters =(1:numSubdivisions-1)/numSubdivisions; 
            PointsOnLine = (1-parameters')*0.5*(PointsIclass(ipoint,:)+ PointsJclass(ipoint,:)) + ...
                0.5*parameters'*(PointsIclass(ipoint+1,:) + PointsJclass(ipoint+1,:));

            if (size(StartPointsLeft, 1) <iidxNextSegment-1)
                StartPointsLeft = [StartPointsLeft; zeros(2*nTriplets-1,ndim)]; %#ok<AGROW>
                StartPointsRight = [StartPointsRight; zeros(2*nTriplets-1,ndim)]; %#ok<AGROW>
            end
            StartPointsLeft(IdxSegment(ipoint):iidxNextSegment-1,:) = PointsOnLine;

            AuxVec(1:ndim) = PointsIclass(ipoint+1,:) - PointsIclass(ipoint,:);
            
            % compute unit normal vector
            NormalVec = cross(AuxVec, NormalVecOfPlane);
            NormalVec = NormalVec/norm(NormalVec(1:ndim));

            % for 2 dimensions: use refined method of computing starting
            % values based upon estimated curvature
            if ndim == 2
                if ipoint > 2 && ipoint < nTriplets-2
                    curv = estimateCurvature(0.5*(PointsIclass(ipoint-2:ipoint+2,:) + ...
                                                 PointsJclass(ipoint-2:ipoint+2,:)), 3, abstolBisection);
            
                % estimating the curvature is most reliable for the
                % midpoint of the subset, therefore, we enlarge its value
                % by an additional safety factor here.
                elseif ipoint <= 2
                    if (nTriplets > 2)
                        curv = safetyFactor*estimateCurvature(0.5*(PointsIclass(1:min(nTriplets,5),:) + ...
                                                                  PointsJclass(1:min(nTriplets,5),:)), ...
                                                             2, abstolBisection);
                    else
                        curv = 0;
                    end
                % in this case, ipoint >= 3, such that we can compute the
                % curvature in any case    
                elseif ipoint >= nTriplets-2
                    curv = safetyFactor*estimateCurvature(0.5*(PointsIclass(max(1, nTriplets-4):nTriplets,:) + ...
                                                              PointsJclass(max(1, nTriplets-4):nTriplets,:)), ...
                                                         4, abstolBisection);
                end
            else
                if ipoint > 2 && ipoint < nTriplets-2
                    curv = estimateCurvature(0.5*(PointsIclass(ipoint-2:ipoint+2,~NormalVecOfPlane) + ...
                                                 PointsJclass(ipoint-2:ipoint+2,~NormalVecOfPlane)), ...
                                            3, abstolBisection);
                elseif ipoint <= 2
                    if (nTriplets > 2)
                        curv = safetyFactor*estimateCurvature(0.5*(PointsIclass(1:min(nTriplets,5),~NormalVecOfPlane) + ...
                                                                  PointsJclass(1:min(nTriplets,5),~NormalVecOfPlane)), ...
                                                             2, abstolBisection);
                    else
                        curv = 0;
                    end
                elseif ipoint >= nTriplets-2
                    curv = safetyFactor*estimateCurvature(0.5*(PointsIclass(max(1,nTriplets-4):nTriplets,~NormalVecOfPlane) + ...
                                                              PointsJclass(max(1,nTriplets-4):nTriplets,~NormalVecOfPlane)), ...
                                                         4, abstolBisection);
                end
            
            end
            
            % Some elaborate computation shows that the maximal deviation
            % from a straight line between two points is
            % 0.25*curv*dist^2 + 1/16*curv^3.*dist^4 + O(h^6) for a planar
            % curve with curvature curv. However, our points are known up
            % to abstolBisection only. Therefore, we should at least
            % consider this deviation from a straight line. If the
            % estimation of the curvature is way too large, we fall back to
            % the old heuristics epsSafemax*dist as deviation.
            errorInd = 0.25*curv*distNext(ipoint)^2 + 1/16*curv^3.*distNext(ipoint)^4;
            errorInd = min(epsSafemax*distNext(ipoint), ...
                           max(0.95*abstolBisection, errorInd));
            
            PointsAux = PointsOnLine + errorInd*NormalVec(1:ndim);
            
            % guarantee that PointsAux is within the domain
            PointsAux = min(max(ProblemDescr.Xmin + eps, PointsAux), ProblemDescr.Xmax - eps);

            StartPointsRight(IdxSegment(ipoint):iidxNextSegment-1,:) = PointsAux;

            PointsAux = PointsOnLine - errorInd*NormalVec(1:ndim);
            PointsAux = min(max(ProblemDescr.Xmin + eps, PointsAux), ProblemDescr.Xmax - eps);
            StartPointsLeft(IdxSegment(ipoint):iidxNextSegment-1,:) = PointsAux;
        end
    end
    
    % shorten and create vectors accordingly
    StartPointsLeft = StartPointsLeft(1:iidxNextSegment-1, :);
    StartPointsRight = StartPointsRight(1:iidxNextSegment-1, :);

    % if there are no points to add at all, skip the computation
    if (size(StartPointsLeft,1) == 0)
        success = true;
        return
    end

    AuxArr = computeClassification([StartPointsLeft; StartPointsRight], ProblemDescr);

    ClassLeft = AuxArr(1: iidxNextSegment-1, :);
    ClassRight = AuxArr(iidxNextSegment:size(AuxArr, 1));

    % find points near the fault line, each with a counterpart
    % in the opposite class
    [PointPairs, IdxSucceeded, ClassPointsSucceeded] = startPairs(StartPointsLeft, ...
        StartPointsRight, ClassLeft, ClassRight, ClassVals(iclass), ClassVals(jclass), ProblemDescr);

    % compute points near the fault line by bisection
    [PointsLeft, PointsRight, Finished] = ...
        computeSurfacePoints(PointPairs(IdxSucceeded,1:ndim), ...
                             PointPairs(IdxSucceeded, ndim+1:2*ndim), ...
                             ClassPointsSucceeded(IdxSucceeded,1), ...
                             ClassPointsSucceeded(IdxSucceeded,2), ...
                             ProblemDescr, FaultApproxParams);

    % Compute array with the indices of the points where computeSurfacePoints
    % suceeded in the complete array of points
    aux = 1:size(IdxSucceeded,1);
    IdxSucceeded = aux(IdxSucceeded);
    
    %avoid to add almost duplicate points
    for i = 1: size(Finished,1)

        if Finished(i)

            % compute distance to all points on the fault line
            DistVec = PointsIclass - PointsLeft(i,:);

            DistVec = DistVec.^2;
            for idim = 2: ndim
                DistVec(:,1) = DistVec(:,1) + DistVec(:,idim);
            end
            DistVec = DistVec(:,1);

            if (sqrt(min(DistVec)) > FaultApproxParams.maxDistForSurfacePoints*FaultApproxParams.minDistFactor)

                nTriplets = nTriplets + 1;

                tripletsAdded = true;

                if (ClassPointsSucceeded(IdxSucceeded(i),1) == iclass)
                    PointsIclass(nTriplets,:) = PointsLeft(i,:);
                    PointsJclass(nTriplets,:) = PointsRight(i,:);      
                else
                    PointsIclass(nTriplets,:) = PointsRight(i,:);
                    PointsJclass(nTriplets,:) = PointsLeft(i,:);      
                end
            end
        end
    end
    success = true;
end