% This function realizes (sub)algorithm adapt in two dimensions. It locally
% refines and coarsens a set of triplets representing the icomp-th
% component of the fault \Gamma_{iclass,jclass} according to the estimated
% error of a polygonal approximation of the fault based on the given set of
% triplets. Hereby, we always rely on the "medium" points of the triplets.
%
% Input:
% - PointSetsSurface: structure of arrays containing the point sets
%   which represent the faults.
% - NumPointsSurf: structure of arrays containing the number of points per
%   set
% - iclass, jclass: class indices
% - icomp: index of the component of \Gamma_{iclass,jclass}
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
% - ClassVals: Array containing the class values. These are not
%   necessarily the class indices. Imagine that f(\Omega) = {1,2,5}. Then,
%   ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
%   Size: nclasses
%
% Output:
% - PointSetsSurface: Structure of arrays containing the point sets
%   which represent the faults.
% - NumPointsSurf: Structure of arrays containing the number of points per
%   set

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointSetsSurface, NumPointsSurf] = adapt2D(PointSetsSurface, ...
    NumPointsSurf, iclass, jclass, icomp, ...
    ProblemDescr, FaultApproxParams, ClassVals)
    
    % aux array for adaptive insertion of triplets
    Iidxtry = [0,1,-1,2,-2];

    % maximal admissible angle between subsequent line segments
    maxAdmissibleAngle = 3.1;
        
    % maximal number of adaptive refinement/coarsening steps
    maxiterAdapt = FaultApproxParams.maxiterAdapt;

    % for safeguarding the search for initial starting pairs for bisection
    epsSafemax = FaultApproxParams.alpha;
    
    % dimension
    ndim = size(PointSetsSurface{iclass,jclass}{icomp}, 2);

    for iiterAdapt = 1: maxiterAdapt
        numPoints = NumPointsSurf{iclass, jclass}(icomp);

        % Euclidean length of each line segment. The line segments are
        % created by connecting the "medium" points of subsequent triplets.
        SegLength = vecnorm((PointSetsSurface{iclass, jclass}{icomp}(1:numPoints-1,:) - ....
                            PointSetsSurface{iclass, jclass}{icomp}(2:numPoints,:))');
        SegLength = SegLength';
        
        if (numPoints > 2)
            % Compute the curvature and new starting points for creating
            % starting pairs in case of refinement.
            [curvature, NewPoints] = CurvAndNewPoints(0.5*(PointSetsSurface{iclass, jclass}{icomp} + ...
                                                           PointSetsSurface{jclass, iclass}{icomp}), ...
                                                      numPoints, FaultApproxParams.abstolBisection);
        else
            break
        end
        
        % We approximate the maximal curvature of a line segment with the
        % maximum of the ones at its start- and end point.
        curvAux = max(curvature(1:end-1), curvature(2:end));

        % The error indicator applies to the approximation with line
        % segments. The approximation with RBFs used for finding starting
        % values for bisection is far better.
        % As heuristic, we just take the higher order error contribution as
        % estimation of the deviation of the true curve from the RBF
        % approximation.
        errorIndHOT = 1/16*curvAux.^3.*SegLength.^4;
        errorInd = 0.25*curvAux.*SegLength.^2 + errorIndHOT;

        % As the position of the "medium" points is accurate up to
        % abstolBisection only, it makes sense to further refine
        % only line segments sufficiently longer than this.
        % The factor 3 is due to Pythagoras and the fact that the distance
        % between two approximations of the same point on the fault line is
        % at most 2*abstolBisection.
        SegsToRefine = (errorInd > FaultApproxParams.errMax) & ...
                        SegLength > 3.0* FaultApproxParams.abstolBisection;

        if (any(SegsToRefine))
            % We aim at inserting new triplets in the approximate middle
            % of selected line segments.

            % Hack for using circular shift (if the last segment is to be
            % refined, this last index is shifted at the beginning of the
            % vector.) We avoid this by adding false at the end. Note that
            % the number of points exceeds the number of line segments by
            % one.
            SegsToRefine = [SegsToRefine; false];
            if numPoints > 3
                MidPoints = NewPoints(SegsToRefine,:);
            else
                MidPoints = 0.5*(PointSetsSurface{iclass, jclass}{icomp}(circshift(SegsToRefine,1),:) + ...
                                 PointSetsSurface{iclass, jclass}{icomp}(SegsToRefine,:));
            end
            normals = (PointSetsSurface{iclass, jclass}{icomp}(circshift(SegsToRefine,1),:) - ...
                       PointSetsSurface{iclass, jclass}{icomp}(SegsToRefine,:));
            normals(:,2) = -normals(:,2);
            auxVec = normals(:,2);
            normals(:,2) = normals(:,1);
            normals(:,1) = auxVec;
            normNormals = sqrt(normals(:,1).*normals(:,1) + normals(:,2).*normals(:,2));
            normals = normals./normNormals;

            % We want to find triplets in the approximate middle of two
            % consecutive triplets on the fault line by bisection. For
            % this, we need starting values. We have a heuristic
            % estimation for how much a point on the true fault line may
            % deviate from its RBF approximation (this is what errIndHOT
            % estimates). Therefore, we can use this information for
            % getting reasonable starting values for bisection. However, we
            % need to safeguard this guess by epsSafemax*dist (if the
            % curvature is largely overestimated) and the tolerance for
            % bisection, as points are on any fault only up to this
            % tolerance.
            % If this heuristic fails, startPairs still finds valid
            % starting pairs in almost all cases, but at the cost of
            % additional function evaluations.
            alpha = min(epsSafemax*SegLength(SegsToRefine), ...
                        max(errorIndHOT(SegsToRefine), ...
                            0.95*FaultApproxParams.abstolBisection)) ;
            PointsRight = MidPoints + normals.*alpha;
            PointsRight = min(max(ProblemDescr.Xmin + eps, PointsRight), ProblemDescr.Xmax - eps);

            PointsLeft = MidPoints - normals.*alpha;
            PointsLeft = min(max(ProblemDescr.Xmin + eps, PointsLeft), ...
                             ProblemDescr.Xmax - eps);

            IidxToRefine = find(SegsToRefine == true);
            AuxArr = computeClassification([PointsRight; PointsLeft], ...
                                           ProblemDescr);

            ClassPointsRight = AuxArr(1: size(PointsRight,1), :);
            ClassPointsLeft = AuxArr(size(PointsLeft,1)+1:size(AuxArr, 1));

            [PointPairs, IdxSucceeded, ClassPointsSucceeded] = startPairs(PointsRight, ...
                PointsLeft, ClassPointsRight, ClassPointsLeft, ...
                ClassVals(iclass), ClassVals(jclass), ProblemDescr);

            [PointsRight, PointsLeft, Finished] = tripletsByBisection(PointPairs(IdxSucceeded,1:ndim), ...
            PointPairs(IdxSucceeded, 3:2*ndim), ClassPointsSucceeded(IdxSucceeded,1), ...
            ClassPointsSucceeded(IdxSucceeded,2), ProblemDescr, FaultApproxParams);

            aux = 1:size(IdxSucceeded,1);
            IdxSucceeded = aux(IdxSucceeded);

            tripletsAdded = 1;
            for i = 1: size(IdxSucceeded,2)
                binsertionSuccessful = false;
                if (Finished(i))
                    % insert new triplets
                    if (ClassPointsSucceeded(i,1) == ClassVals(iclass))
                        NewPointIclass = PointsLeft(i,:);
                        NewPointJclass = PointsRight(i,:);
                    elseif (ClassPointsSucceeded(i,1) == ClassVals(jclass))
                        NewPointIclass = PointsRight(i,:);
                        NewPointJclass = PointsLeft(i,:);
                    
                    % It may happen that a new valid triplet has been
                    % found, which however belongs to another fault line.
                    % Then, we skip this triplet and move on to the next
                    % one. This part of the implementation should be
                    % improved.
                    else
                        continue
                    end
                    
                    % It is not given that if we intend to insert a new
                    % triplet in between two consecutive triplets i and
                    % i+1, that the new triplet found is indeed in between
                    % these consecutive triplets.
                    % Therefore, we test this. If this is not the case, we
                    % test some of the neighbouring line segments for
                    % insertion if they exist. If this still fails, we give
                    % up.
                    for itest = 1:size(Iidxtry,2)
                        ishift = Iidxtry(itest);

                        iidxStart = IidxToRefine(IdxSucceeded(i))+tripletsAdded+ishift-1;
                        
                        % if the corresponding line segment exists
                        if NumPointsSurf{iclass, jclass}(icomp) >= iidxStart+1 && iidxStart > 0
                            PointsTest = [PointSetsSurface{iclass, jclass}{icomp}(1: iidxStart,:); ...
                                          NewPointIclass; ...
                                          PointSetsSurface{iclass, jclass}{icomp}(iidxStart+1:NumPointsSurf{iclass, jclass}(icomp),:)];

                            % We want to insert a triplet "somewhere" in
                            % the middle between two subsequent triplets
                            % with indices iidxStart and iidxStart+1.
                            % However, we need to ensure that the new
                            % triplet is indeed somewhere in the middle.
                            % For doing so, we compare the distance of the
                            % new triplet to the given ones relative to
                            % their distance.
                            % If minDist is smaller than some threshold,
                            % then the new triplet is very close to one of
                            % the existing triplets.
                            minDist = min(norm(PointSetsSurface{iclass, jclass}{icomp}(iidxStart,:) - NewPointIclass,2), ...
                                      norm(PointSetsSurface{iclass, jclass}{icomp}(iidxStart+1,:) - NewPointIclass,2));

                            minDist = minDist/norm(PointSetsSurface{iclass, jclass}{icomp}(iidxStart,:) - ...
                                                   PointSetsSurface{iclass, jclass}{icomp}(iidxStart+1,:));
                                  
                            % If minDist is not in the range given here,
                            % the new triplet is not in the middle at all,
                            % such that we do not insert it.
                            % The range is usually violated in case of
                            % wrong sorting only.
                            if minDist > 0.2 && minDist < 1.0
                                maxAngle = computeMaximalAngle(PointsTest(max(1,iidxStart-2): min(iidxStart+3,size(PointsTest,1)),:));

                                noIntersection = ~selfIntersection(PointsTest);

                                if noIntersection && (maxAngle < maxAdmissibleAngle)
                                    PointSetsSurface{iclass, jclass}{icomp} = PointsTest;
                                    PointSetsSurface{jclass, iclass}{icomp} = ...
                                        [PointSetsSurface{jclass, iclass}{icomp}(1: iidxStart,:); ...
                                         NewPointJclass; ...
                                         PointSetsSurface{jclass, iclass}{icomp}(iidxStart+1:NumPointsSurf{jclass, iclass}(icomp),:)];

                                    % If we have found a suitable line
                                    % segment, we can stop.
                                    binsertionSuccessful = true;
                                    break
                                end
                            end
                        end
                    end

                    if binsertionSuccessful 
                        NumPointsSurf{iclass, jclass}(icomp) = NumPointsSurf{iclass, jclass}(icomp) + 1;
                        NumPointsSurf{jclass, iclass}(icomp) = NumPointsSurf{jclass, iclass}(icomp) + 1;
                        tripletsAdded = tripletsAdded + 1;
                    else
                        warning(['adding triplets in adaptive refinement failed for a fault line between classes ', ...
                                 int2str(iclass), ' and ', int2str(jclass)])
                    end
                end
            end

            numPoints = NumPointsSurf{iclass, jclass}(icomp);
            
            SegLength = vecnorm((PointSetsSurface{iclass, jclass}{icomp}(1:NumPointsSurf{iclass, jclass}(icomp)-1,:) - ....
                                 PointSetsSurface{iclass, jclass}{icomp}(2:NumPointsSurf{iclass, jclass}(icomp),:))');
            SegLength = SegLength';
            
            if (numPoints > 2)
                % curvature on any triplet
                curvature = CurvAndNewPoints(0.5*(PointSetsSurface{iclass, jclass}{icomp} + ...
                                                  PointSetsSurface{jclass, iclass}{icomp}), ...
                                             numPoints, FaultApproxParams.abstolBisection);

                % We approximate the maximal curvature of a line segment
                % with the maximum of the ones at its start- and end point.
                curvAux = max(curvature(1:end-1), curvature(2:end));

            else
                break
            end

            % If points have been added, recompute the error indicator per
            % line segment.
            errorInd = 0.25*curvAux.*SegLength.^2 + 1/16*curvAux.^3.*SegLength.^4;
        end

        % --- coarsening step ---

        % Delete a triplet, if the estimated error of both line segments
        % the triplet belongs to is smaller than errMin.
        TipletsToRemove = errorInd < FaultApproxParams.errMin;

        % Heuristic: Never delete the first or the last triplet.
        TipletsToRemove = and(TipletsToRemove(1:end-1), TipletsToRemove(2:end));
        TipletsToRemove = [0; TipletsToRemove; 0];

        if (~any(any(SegsToRefine)) && ~any(TipletsToRemove))
            break
        end

        % Heuristic: Never delete consecutive triplets.
        for i = 2:size(TipletsToRemove,1)-1
            TipletsToRemove(i) = TipletsToRemove(i) & ~TipletsToRemove(i-1); 
        end

        PointSetsSurface{iclass,jclass}{icomp} = PointSetsSurface{iclass,jclass}{icomp}(~TipletsToRemove,:);
        PointSetsSurface{jclass,iclass}{icomp} = PointSetsSurface{jclass,iclass}{icomp}(~TipletsToRemove,:);
        NumPointsSurf{iclass, jclass}(icomp) = size(PointSetsSurface{jclass,iclass}{icomp},1);
        NumPointsSurf{jclass, iclass}(icomp) = NumPointsSurf{iclass, jclass}(icomp);
    end
    

    function [curvature, StartPoints] = CurvAndNewPoints(PointSet, numPoints, abstolBisection)
        curvature = zeros(numPoints, 1);
        StartPoints = zeros(numPoints-1, 2);

        % For all "inner" points: consider five points for computing the
        % normal.
        % This applies for all components containing at least five points.
        for ipoint = 3: numPoints-2
            PointsLocal = PointSet(ipoint-2:ipoint+2,:);
            [curvature(ipoint), startval] = estimateCurvature(PointsLocal, 3, abstolBisection);
            StartPoints(ipoint-1:ipoint,:) = StartPoints(ipoint-1:ipoint,:) + startval;
        end                            

        StartPoints(3:end-2,:) = 0.5*StartPoints(3:end-2,:);

        % Consider the first, second, last and second but last points, if
        % feasible.
        if (numPoints >= 4)
            PointsLocal = PointSet(1:min(5,numPoints),:);
            curvature(1) = estimateCurvature(PointsLocal, 1, abstolBisection);
            [curvature(2), startval] = estimateCurvature(PointsLocal, 2, abstolBisection);
            StartPoints(1,:) = startval(1,:);

            PointsLocal = PointSet(end-min(5,numPoints)+1:end,:);
            [curvature(end-1), startval] = estimateCurvature(PointsLocal, min(5,numPoints)-1, abstolBisection);
            StartPoints(end,:) = startval(2,:);
            curvature(end) = estimateCurvature(PointsLocal, min(5,numPoints), abstolBisection);
            
        % Fault component consists of three points only: compute the
        % curvature based upon these three points.
        elseif (numPoints == 3)
            PointsLocal = PointSet;
            curvature(1) = estimateCurvature(PointsLocal, 1, abstolBisection);
            [curvature(2), StartPoints] = estimateCurvature(PointsLocal, 2, abstolBisection);
            curvature(3) = estimateCurvature(PointsLocal, 3, abstolBisection);

        % Two points only: We can not obtain any curvature information, so
        % we just return 0 and as starting value the midpoint.
        else
            curvature = [0;0];
            StartPoints = 0.5*(PointSet(1,:) - PointSet(2,:));
        end
    end
end