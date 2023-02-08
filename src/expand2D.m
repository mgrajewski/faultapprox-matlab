% This function extends a fault line up to its end. We assume that the
% points on the fault line are ordered. There are two modes of
% operation: mode 0 for extending the fault line before its current
% starting point and mode 1 for extending it beyond its current end.
% Depending on mode, we take the first or last three points (if
% existing) of the fault line and fit an interpolating polynomial
% curve. We then create new points on this line such that their
% distance is approximately the distance of the ones already existing.
% These points and the corresponding auxiliary points are used to
% obtain new points on the fault line.
% We repeat this process until we leave the domain or we
% intrude in a subdomain neither belonging to iclass nor jclass. In
% this case, we take the first of these points and discard all others
% of this kind. We adjust its parameter by bisection to obtain the
% maximal parameter where no third class is involved. This point marks
% approximately the end of the fault line.
% For any fault line, we store at which side of the domain the fault
% line leaves the domain if so in the array LeftDomain. This information
% is necessary for constructing the approximating polygons of any
% subdomain in function ReconstructSubdomains. The information is coded
% as follows: 1: left boundary, 2: right boundary, 3: bottom boundary,
% 4: top boundary, 0: fault line end inside the domain
% as follows: 1: bottom boundary, 2: right boundary, 3: top boundary,
% 4: left boundary, 0: fault line end inside the domain
%
% This routine is applied for expanding fault surfaces in 3D as well.
% Here, we take some points from the fault surface and treat them as
% points defining a curve to continue. Again, we need points near the
% curve on both sides of the fault line. Using just a vector normal to
% the curve as in 2D does not work. Therefore, we need additional
% information which is provided by NormalVecPlane.
%
% Inoutput:
% - LeftDomain: description see above
% - npoints: number of points on the fault lines
% - PointsIclass, PointsJclass: the fault points near the current fault
%   line in either class iclass or class jclass
% 
% Input:
% IdxPoints: array of point indices according to
%   their position on the fault line
% - iclass, jclass: class indices
% - icomp: index of the fault component 
% - ClassVals: Array containing the class values. These are not
%   necessarily the class indices. Imagine that f(\Omega) = {1,2,5}. Then,
%   ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
%   Size: nclasses
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
% - mode: mode of operation (0 or 1, see description above)
%
% Output:
% - IdxPointsEnlarged: arrray of point indices according to
%   their position on the fault line, but enlarged by the indices of the
%   points added
% - resort: If the fault line is closed, the start and the end of the
%   fault line in fact overlap. This requires a complete resorting of the
%   points. We do not resort inside this routine but indicate the
%   requirement by setting resort = true.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [IdxPointsEnlarged, LeftDomain, npoints, PointsIclass, ...
          PointsJclass, resort, pointsAdded] = ...
          expand2D(LeftDomain, IdxPoints, iclass, jclass, icomp, ...
                   ClassVals, PointsIclass, PointsJclass, ...
                   NormalVecPlane, ProblemDescr, FaultApproxParams, mode)


    numPointsOnCurve = 4;
    pointsAdded = false;

    % maximal number of expansion cycles (to prevent endless loops)
    nmaxcycles = 500;

    finished = false;
    bdomainLeft = false;

    % distVec contains actual distance between consecutive points
    numPoints = size(PointsIclass, 1);

    if numPoints < 2
        warning('Not enough points provided for expanding the fault line.')
        IdxPointsEnlarged = IdxPoints;
        npoints = numPoints;
        resort = false;
        pointsAdded = false;
        return
    end

    distVecAux = (PointsIclass(IdxPoints(1:numPoints-1),1) - ...
                  PointsIclass(IdxPoints(2:numPoints),1)).^2 + ...
                 (PointsIclass(IdxPoints(1:numPoints-1),2) - ...
                  PointsIclass(IdxPoints(2:numPoints),2)).^2 ;

    distVecAux = sqrt(distVecAux);
    distVec(IdxPoints) = [distVecAux(1); ...
                          0.5*(distVecAux(1:numPoints-2) + distVecAux(2:numPoints-1)); ...
                          distVecAux(numPoints-1)];
    
    % requirement for complete resorting (not necessary by default)
    resort = false;

    % for computing valid point pairs as starting values for bisection
    alpha = FaultApproxParams.alpha;

    % points closer as epstol are regarded as identical
    epstol = FaultApproxParams.eps;
    
    [npoints, ndim] = size(PointsIclass);
        
    % auxiliary vector
    AuxVecLoc = [0,0,0];
    
    % index array for the extended set points on the surface (at the
    % beginning, this is the original set of points)
    IdxPointsEnlarged = IdxPoints;
    
    % number of cycles in expansion
    ncycles = 0;

    while ~finished && ncycles <= nmaxcycles

        ncycles = ncycles + 1;

        % expand at the beginning
        if (mode == 0)
            IdxSet = min(npoints, numPointsOnCurve):-1:1;

        % expand at the end
        elseif (mode == 1)
            IdxSet = max(1, npoints - (numPointsOnCurve-1)): npoints;
        end
        numPointsToInterpolate = size(IdxSet,2);
        PointsOnFault = 0.5*(PointsIclass(IdxPointsEnlarged(IdxSet), :) + ...
                             PointsJclass(IdxPointsEnlarged(IdxSet), :));

        % extrapolate current points
        [NewPointOnFault, avg_dist, extraParam, DataPars, Q, xmean, coeffs] = ...
            pointsOnCurveByExtrapolation(PointsOnFault, 1, FaultApproxParams);            


        finishedSurfPoint = false;
        itry = 0;
        while itry <= 3 && ~finishedSurfPoint
            itry = itry + 1;

            extraParam = 1/2^(itry-1)*(extraParam - ...
                                       DataPars(numPointsToInterpolate)) + ...
                                       DataPars(numPointsToInterpolate);

            NewPointOnFault = evalPol(extraParam, coeffs, Q, xmean);
  
            % distance of the new point relative to the average distance of
            % the existing points
            relDist = norm(NewPointOnFault - PointsOnFault(numPointsToInterpolate,:) ,2)/avg_dist;
            
            % test, if the extrapolated point is outside the domain or not
            isOutsideDomain = any(NewPointOnFault < ProblemDescr.Xmin | ...
                                  NewPointOnFault > ProblemDescr.Xmax);
    
            % compute normal vector
            AuxVecLoc(1:ndim) = NewPointOnFault - PointsOnFault(numPointsToInterpolate,:);
            NormalVecLoc = cross(AuxVecLoc, NormalVecPlane);
            NormalVecLoc = NormalVecLoc(1:ndim)/norm(NormalVecLoc, 2);          
            % this is purely heuristical
            stepSize = max(0.95*FaultApproxParams.abstolBisection, ...
                           0.5* min(alpha*avg_dist, ...
                                    3*FaultApproxParams.abstolBisection*(relDist.^2+1)));
                    
            PointLeft = NewPointOnFault + stepSize*NormalVecLoc;
            PointLeft = min(max(ProblemDescr.Xmin + epstol, PointLeft), ...
                            ProblemDescr.Xmax - epstol);
    
            PointRight = NewPointOnFault - stepSize*NormalVecLoc;
            PointRight = min(max(ProblemDescr.Xmin + epstol, PointRight), ...
                             ProblemDescr.Xmax - epstol);
            
            if (~isOutsideDomain)
                AuxArr = computeClassification([PointLeft; PointRight], ProblemDescr);
                
                classPointLeft = AuxArr(1);
                classPointRight = AuxArr(2);
    
                isOutsideSub = (classPointLeft ~=ClassVals(iclass) & ...
                                classPointLeft ~= ClassVals(jclass)) | ...
                               (classPointRight ~=ClassVals(iclass) & ...
                                classPointRight ~= ClassVals(jclass));        
            else
                classPointLeft = -1;
                classPointRight = -1;
    
                isOutsideSub = false;
            end
 
            % If there are some points on the curve which are neither in
            % iclass nor in jclass, we have reached the boundary of the
            % domain or the subdomain and have therefore reached the last
            % cycle of the expansion loop.
            % As points outside the valid domain have class -1, this
            % includes the case the domain has been left.
            % Note that if the "central" points from which the left and the
            % right points was generated is outside the domain, we
            % automatically assume that this holds for both the left and the
            % right points. Concentrating here on just PointsLeft therefore
            % does not harm.
            if isOutsideDomain
                finished = true;
                bdomainLeft = true;
    
                % approximate the intersection of the extrapolated line and
                % the domain boundary.
                % We know that DataPars(numPointsToInterpolate) belongs to
                % a point inside the domain. Therefore, it is a lower bound
                % for the parameter of the intersection. On the other hand,
                % point corresponding to extraParam does not belong to the
                % domain. Therefore, this is an upper bound for the
                % parameter.
                [xnew, iedge] = extraNewPointOnDomainBdry(DataPars(numPointsToInterpolate), ...
                    extraParam, PointsOnFault(end,:), NewPointOnFault, coeffs, Q, xmean, ...
                    ProblemDescr, FaultApproxParams);
                
                % updateVec is not necessarily the normal vector of the
                % extrapolated fault line in xnew, so we call it updateVec.
                updateVec = cross(NormalVecPlane, AuxVecLoc);
                                
                % from xnew, we need to create a pair of points as
                % starting point for bisection. We need to ensure that the
                % line segment between these two points is parallel to the
                % appropriate domain boundary.
                switch iedge
                    
                    % x2 is minimal (bottom boundary/front surface of cube)
                    case 1
                    updateVec(2) = 0;
    
                    % x1 is maximal (right boundary/right surface of cube)
                    case 2
                    updateVec(1) = 0;
    
                    % x2 is maximal (top boundary/back surface of cube)
                    case 3
                    updateVec(2) = 0;
    
                    % x1 is minimal (left boundary/left surface of cube)
                    case 4
                    updateVec(1) = 0;
    
                    % x3 is minimal (bottom surface of cube)
                    case 5
                    updateVec(3) = 0;
    
                    % x3 is maximal (top surface of cube)
                    case 6
                    updateVec(3) = 0;
                end
                
                relDist = norm(xnew - PointsOnFault(numPointsToInterpolate,:))/avg_dist;
                updateVec = updateVec(1:ndim);
                
                %heuristic control of step length
                updateVec = updateVec/norm(updateVec)* ...
                    max(0.95*FaultApproxParams.abstolBisection, ...
                        3*FaultApproxParams.abstolBisection*(relDist^2+1));
                
                xnewLeft = xnew + updateVec;
                xnewRight = xnew - updateVec;
                
                % ensure that the new points near the boundary are really
                % inside the domain
                xnewLeft = max(min(xnewLeft, ProblemDescr.Xmax - epstol), ...
                               ProblemDescr.Xmin + epstol);
                xnewRight = max(min(xnewRight, ProblemDescr.Xmax - epstol), ...
                                ProblemDescr.Xmin + epstol);
    
                % no easier way to do this, as computeClassification
                % outputs just one array, which can not be assigned to
                % classPointLeft and classPointRight directly.
                classNew = computeClassification([xnewLeft; xnewRight], ProblemDescr);
                classPointLeft = classNew(1);
                classPointRight = classNew(2);
                LeftDomain(iclass, jclass, icomp) = iedge;
    
                % find starting values for bisection
                [PointPairs, IdxOk, ClassPointsOk] = ...
                    startPairs(xnewLeft, xnewRight, ...
                               classPointLeft, classPointRight, ...
                               ClassVals(iclass), ClassVals(jclass), ...
                               ProblemDescr);
                
                if IdxOk
                    [PointLeft, PointRight, finishedSurfPoint] = ...
                        computeSingleSurfacePoint(PointPairs(1:ndim), ...
                        PointPairs(ndim+1:2*ndim), ClassPointsOk(1), ...
                        ClassPointsOk(2), ProblemDescr, FaultApproxParams);
                end
    
                % There is no indication that the fault line has left the
                % domain. Therefore, we can continue to expand unless it
                % ends inside the domain. 
            elseif isOutsideSub
                finished = true;
                bdomainLeft = false;
    
                % Did we expand the fault line beyond its end? If so,
                % either PointLeft or PointRight belongs neither to iclass
                % nor to jclass.
                % In this case, take last known point pair in the subdomain
                % and the first known outside, discard all others and find
                % the end of the extrapolated fault line by bisection.
                % However, this is a heuristic approach, as it may happen
                % that e.g. due to insufficient distance to the
                % extrapolation curve, both points belong to the same
                % class.
                %                      II
                %        -------------------------->     x PointsLeft     
                % ______x________x________x_____         x extrapolated point
                %                               \--__    x PointsRight
                %              I                 \   ---__    
                %                                 \       ---__
                %                                  \    III    ---__
                %
                % Then, we do not detect that we extrapolated the fault
                % line beyond its end. We will fail to find a valid
                % starting pair such that expand stops and no point
                % is added, even if the ed of the fault line is yet not
                % reached.
    
                % A fault line ends inside the domain by intersecting a
                % subdomain with a third class.
                % We aim at finding the intersection of the extrapolated
                % fault line with that third subdomain using bisection.
    
                % A lower bound tmin for the parameter of the intersection
                % is the parameter of the last known point inside the
                % subdomain.
                tmin = DataPars(numPointsToInterpolate);
                
                % Correspondingly, the parameter of the first point on
                % the extrapolated fault line which is known to be
                % outside the subdomain is an upper bound for the
                % parameter of the intersection.
                tmax = extraParam(1);
                xmin = evalPol(tmin, coeffs, Q, xmean);
                xmax = evalPol(tmax, coeffs, Q, xmean);
                dist = norm(xmin - xmax);
    
                xminAuxAll = [];
                xminAll = [];
                ClassMinAll = [];
    
                iiter = 1;
                while iiter < 10 && dist > FaultApproxParams.abstolBisection
    
                    tnew = 0.5*(tmin+ tmax);
                    xnew = evalPol(tnew, coeffs, Q, xmean);
    
                    AuxVecLoc(1:ndim) = xnew - xmin;
                    NormalVecLoc = cross(AuxVecLoc, NormalVecPlane);
                    NormalVecLoc = NormalVecLoc(1:ndim)/norm(NormalVecLoc, 2);
                    relDist = norm(PointsOnFault(numPointsToInterpolate,:) - xnew,2)/avg_dist;
    
                    % this is purely heuristical
                    stepSizes = max(0.95*FaultApproxParams.abstolBisection, ...
                                    0.5* min(alpha*avg_dist, ...
                                             3*FaultApproxParams.abstolBisection*(relDist.^2+1)));
                    xnewLeft = xnew + stepSizes*NormalVecLoc;
                    xnewRight = xnew - stepSizes*NormalVecLoc;
    
                    xnewLeft = min(max(ProblemDescr.Xmin + epstol, xnewLeft), ...
                                   ProblemDescr.Xmax - epstol);
                    xnewRight = min(max(ProblemDescr.Xmin + epstol, xnewRight), ...
                                    ProblemDescr.Xmax - epstol);
    
                    classNew = computeClassification([xnewLeft; xnewRight], ProblemDescr);
    
                    if (any(classNew ~= ClassVals(iclass) & classNew ~= ClassVals(jclass)))
                        tmax = tnew;
                        xmax = xnew;
                    else
                        tmin = tnew;
                        xmin = xnew;
                        xminAuxAll = [xminAuxAll; [xnewLeft, xnewRight]];
                        xminAll = [xminAll; xmin];
                        ClassMinAll = [ClassMinAll; classNew'];
                    end
    
                    dist = norm(xmax - xmin);
                    iiter = iiter + 1;
                end
    
                iiter = 0;
                finishedSurfPoint = false;
                while iiter < size(xminAll,1) && ~finishedSurfPoint
                    iiter = iiter + 1;
                    xnew = xminAll(end-iiter+1,:);
    
                    % We need a point which is inside the current subdomain,
                    % and we know that xmin is inside the domain. However,
                    % if tmin is (almost) the parameter value of the last
                    % known point, xmin is in fact the last known point on
                    % the fault. It does not make sense to add this point
                    % again.
                    if (norm(xnew - PointsOnFault(numPointsToInterpolate,:), 2))> FaultApproxParams.abstolBisection
        
                        xnewLeft = xminAuxAll(end, 1:ndim);
                        xnewRight = xminAuxAll(end, ndim+1:2*ndim);
                
                        classNew = zeros(1,2);
                        classNew(1) = computeClassification(xnew, ProblemDescr);
                        
                        xnewSuitable = classNew(1) == ClassVals(iclass) || classNew(1) == ClassVals(jclass);
    
    
                        if xnewSuitable
                            % xnew belongs to a valid point pair
                            if (classNew(1) == ClassVals(iclass) && ClassMinAll(end-iiter+1, 1) == ClassVals(jclass))
                                xnewLeft = xnew;
                                classPointLeft = classNew(1);
                                xnewRight = xminAuxAll(end-iiter+1, 1:ndim);
                                classPointRight = ClassMinAll(end-iiter+1,1);
                            elseif (classNew(1) == ClassVals(jclass) && ClassMinAll(end-iiter+1, 1) == ClassVals(iclass))
                                xnewLeft = xnew;
                                classPointLeft = classNew(1);
                                xnewRight = xminAuxAll(end-iiter+1, 1:ndim);
                                classPointRight = ClassMinAll(end-iiter+1,1);
                            elseif (classNew(1) == ClassVals(jclass) && ClassMinAll(end-iiter+1, 2) == ClassVals(iclass))
                                xnewLeft = xnew;
                                classPointLeft = classNew(1);
                                xnewRight = xminAuxAll(end-iiter+1, ndim+1:2*ndim);
                                classPointRight = ClassMinAll(end-iiter+1,2);
                            elseif (classNew(1) == ClassVals(iclass) && ClassMinAll(end-iiter+1, 2) == ClassVals(jclass))
                                xnewLeft = xnew;
                                classPointLeft = classNew(1);
                                xnewRight = xminAuxAll(end-iiter+1, ndim+1:2*ndim);
                                classPointRight = ClassMinAll(end-iiter+1,2);
                            end
            
                            bdomainLeft = true;
                            [PointLeft, PointRight, finishedSurfPoint] = computeSingleSurfacePoint(xnewLeft, ...
                            xnewRight, classPointLeft, ...
                            classPointRight, ProblemDescr, FaultApproxParams);
                            ClassPointsOk = [classPointLeft, classPointRight];
                        end    
                    % attach boundary point
                    end
                end
                LeftDomain(iclass, jclass, icomp) = 0;
    
                
            % expansion continues
            else
    
                % find starting values for bisection
                [PointPairs, IdxOk, ClassPointsOk] = startPairs(PointLeft, ...
                    PointRight, classPointLeft, classPointRight, ClassVals(iclass), ClassVals(jclass), ProblemDescr);
    
                % If no point has been found at all, stop. If we have no indication
                % that the fault line left the domain, we must assume that it ended
                % inside.
                if IdxOk
                    [PointLeft, PointRight, finishedSurfPoint] = computeSingleSurfacePoint(PointPairs(1:ndim), ...
                    PointPairs(ndim+1:2*ndim), ClassPointsOk(1), ...
                    ClassPointsOk(2), ProblemDescr, FaultApproxParams);
                else
                    LeftDomain(iclass, jclass, icomp) = 0;
                    finished = true;
                end
            end
    
            % add new point if appropriate
            if finishedSurfPoint
    
                % Compute the distance to currently the last point on the
                % fault line If this distance is extremely small, replace
                % this point by the new one.
                LastCurrentPoint = PointsOnFault(end,:);
                if ClassPointsOk(1) ==ClassVals(iclass)
                    distToNewPoint = norm(LastCurrentPoint - PointLeft,2);
                else
                    distToNewPoint = norm(LastCurrentPoint - PointRight,2);
                end
    
                % the new point is sufficiently far away from the last
                % known one: add it
                if (distToNewPoint > avg_dist*FaultApproxParams.minDistFactor)
                    npoints = npoints + 1;
                    iidxAdd = npoints;
    
                    if (mode == 0)
                        IdxPointsEnlarged = [npoints, IdxPointsEnlarged];
                    elseif (mode == 1)
                        IdxPointsEnlarged = [IdxPointsEnlarged, npoints];
                    end            
                
                % the new point is almost the last known one: replace the
                % last known one by the new one and stop expanding
                else
                    iidxAdd = IdxPointsEnlarged(IdxSet(end));
                    finished = true;
                end
                
                if (ClassPointsOk(1) ==ClassVals(iclass))
                    PointsIclass(iidxAdd,:) = PointLeft;
                    PointsJclass(iidxAdd,:) = PointRight;
                else
                    PointsIclass(iidxAdd,:) = PointRight;
                    PointsJclass(iidxAdd,:) = PointLeft;
                end
    
                % If the fault line appears to start or end inside the
                % domain, it may be closed. We test this here.
                if LeftDomain(iclass, jclass, icomp) <= 0
                    closedComp = testForClosedComp(PointsIclass, ...
                                                   IdxPointsEnlarged, ...
                                                   distVec, npoints, mode);      
        
                    if (closedComp)
                        % If the fault component is closed, it does not 
                        % start or end at the domain boundary. We finish
                        % expanding this fault line and state that it
                        % starts and ends inside the domain.
                        LeftDomain(iclass, jclass, icomp) = 0;
        
                        % The new point is in fact somewhere between
                        % existing points. This requires resorting.
                        resort = true;
                        finished = true;
                    end
                end
    
                pointsAdded = true;
                distVec(iidxAdd) = distToNewPoint;
    
            % no suitable point for adding found: try again. We assume
            % in this case that the fault line starts or ends inside the
            % domain unless found otherwise.
            else
                if (LeftDomain(iclass, jclass, icomp) < 0 && ~bdomainLeft)
                    LeftDomain(iclass, jclass, icomp) = 0;
                end
                %finished = true;
            end      
        end
    end
end