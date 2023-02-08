% This function sorts points on a fault line based on distances. We
% assume at first that the points on the fault belong to one fault line
% only. We furthermore assume that a fault line is continuous. It may
% contain kinks, but at these, the angle must not exceed
% acos(FaultApproxParams.cosAngleMax). If we start with the first point
% on the line, the next one is the nearest to the first, unless it is
% included in the set of already sorted points or if the angle between
% the current line segment and the one with the new point exceeds
% acos(FaultApproxParams.cosAngleMax). Then, we proceed with the
% second-nearest point, and so on.  on. If this fails for all
% knearestNeighbours nearest points, we stop sorting. However, there
% will be leftover points. We restart the algorithm with these
% leftover points. If there are no leftover points any more, we are
% finished.
% Leftover points may occur either if we do not start with first or
% last points on the fault or if the fault consists of several
% components (this may happen if subdomains are not simply connected).
% Finding the correct start or end point is hard, and we use a heuristic
% based on the distance to the domain boundary.
% This way, we end up with a bunch of ordered subsets. Ideally, this
% corresponds to the fact that the fault line consists of several
% components, one subset for each component.
%
% In the last step, we try to combine these subsets based (again) on
% distances, if bforceCombine is set to true.
% We made this last step optional, as it may happen that a fault
% consists of several components. Then it does not make sense to sew
% them together.
%
% Input:
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - IdxPointsSurfOrdered: index rray of the sorted point set
% - sortingSuccessful: true, if the points could be sorted, false,
%   otherwise

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [IdxPointsSurfOrdered, sortingSuccessful]  = ...
    sortPointsOnFaultLine(PointSet, bforceCombine, ProblemDescr, ...
                          FaultApproxParams)
    
    knearestNeighbours = FaultApproxParams.nNearestSorting;
    
    ncomps = 0;
    numPoints = size(PointSet,1);

    sortingSuccessful = false;
    
    % if there are two points or even less, there is no need for sorting
    if (numPoints <= 2)
        IdxPointsSurfOrdered{1} = 1:numPoints;
        sortingSuccessful = true;
        return
    end
    
    % heuristic number of trials aka the maximal number of sorted subsets
    % If there are more than about half of the number of points, we must
    % assume that the sorting did not yield any reasonable result.
    ntrials = max(2, floor(numPoints*0.5));
    
    % number of points per sorted subset
    numPointsPerPart = zeros(1, ntrials);
    ndim = size(PointSet,2);

    % as PointSet will be reduced to the number of still unsorted points,
    % we should store the original point set.
    PointSetOrig = PointSet;
    
    IdxPointsSurfOrdered = cell(1);
    
    % compute distance matrix
    DistMatAux = zeros(numPoints, numPoints, ndim);
    for k = 1: numPoints
        DistMatAux(:,k,:) = PointSet - PointSet(k,:);
    end

    DistMatAux = DistMatAux.^2;
    for idim = 2: ndim
        DistMatAux(:,:,1) = DistMatAux(:,:,1) + DistMatAux(:,:,idim);
    end
    DistMat = DistMatAux(:,:,1).^0.5;
    [~, IidxNearestNeighbour] = sort(DistMat);
    IdxAux = 1:numPoints;
    
    
    for itry = 1:ntrials

        % All what follows makes only sense if there is more than one
        % point. If numPoints is one, then there is one orphaned point,
        % and we make it a new segment.
        % Note that numPoints is reduced during the number of trials.
        if (numPoints > 1)
            % Search the points closest to the beginning or end of the fault
            % line.
            % In general, it is impossible to decide which of the point is
            % closest to the beginning/end without knowing the fault line. As
            % heuristic measure, we search the points closest to boundaries
            % of the domain.
            %
            % auxArr contains the distances to all 4 (in 2D) or 6 (in 3D)
            % lines/planes the domain boundary consists of  (even indices
            % of dim for minima, odd indices of dim for maxima)
            auxArr = zeros(2*ndim, numPoints);
            IdxPointsSurfOrdered{itry} = zeros(1,numPoints);
            for idim = 1:ndim
                auxArr(2*idim-1,:) = abs(PointSet(:,idim)' - ProblemDescr.Xmin(idim));
                auxArr(2*idim  ,:) = abs(PointSet(:,idim)' - ProblemDescr.Xmax(idim));
            end
        
        
            % minimal distance to boundary
            [auxArr2, minIndex] = min(auxArr);

            % find all points with minimal distance up to 2*abstolBisection
            Candidates = auxArr2 < min(auxArr2) + 2*FaultApproxParams.abstolBisection;
        
            %      _______________________
            %     |                       |
            %     |                       |
            %     |                       |
            %     |                       |
            %     |                       |
            %     |       |               |
            %     |       x               |
            %     |       |               |
            %     |       x--x-x---x      |
            %     |_______________________|
            % Find out if some points have almost the same distance to the
            % same nearest domain boundary part. In this case, we cannot
            % use this distance for finding a reasonable starting point:
            % For the example above, any of the four bottom points is
            % nearest, depending on numerical inaccuracies.
            % Therefore, we consider these candidates and exclude the
            % minimal distance to the nearest boundary part. Then, we
            % repeat this process.
            %
            % The loop is necessary, as in 3D, points can be along a line
            % parallel to one of the coordinate axes. In this case, we need
            % to exclude 2 planes, as the distance to them of all these
            % points is almost the same and therefore meaningless.
            %
            % It may even happen, that there are several points with almost
            % the same minimal distance to the boundary, but with respect
            % to different lines/planes:
            %      _______________________
            %     |                       |
            %     |                       |
            %     |                       |
            %     |                       |
            %     |                       |
            %     |  |                    |
            %     |  x                    |
            %     |  |                    |
            %     |  x--x-x---x           |
            %     |_______________________|
            %
            % This case is however very special and hard to handle:
            % Considering distances, we would end up with the corner point,
            % which is not the right one. Moreover, sorting works in most
            % cases even for starting points not being the first or last
            % point, so we exclude this special case by 
            % size(unique(minIndex(Candidates)),2) == 1. That means that
            % for the above case, the starting point is more or less
            % arbitrary. However, TestCaseFaultDetection03 shows that even
            % then, sorting works.            
            iskip = 1;
            while (iskip < ndim && size(unique(minIndex(Candidates)),2) == 1 && size(minIndex(Candidates),2) > 1)
                % exclude corresponding coordinate and all non-candidates
                % for starting points
                auxArr = auxArr([1:minIndex(1)-1 minIndex+1:2*ndim-iskip+1], :);

                % setting the distance to some high value is easier
                % than index arithmetics, however, not so efficient
                auxArr(:, ~Candidates) = 1e20;

                % minimal distance to boundary
                [auxArr2, minIndex] = min(auxArr);

                Candidates = auxArr2 < min(auxArr2) + 2*FaultApproxParams.abstolBisection;
                iskip = iskip + 1;
            end
                   
            [~, IdxNearestToBoundary] = sort(auxArr2);
            iidxPoint = IdxNearestToBoundary(1);

            % Ordering points by nearest neighbour search does not work in
            % every case, not even when starting at the beginning or the
            % end of the fault line. However, starting at these points is
            % beneficial, as there are less failure modes. This means that
            % we need a fallback for failed sorting.
            % If the ordering fails for our starting point, we take the
            % second nearest to the boundary and so on.
            IdxPointsSurfOrdered{itry}(1) = iidxPoint;
            ListOfForbiddenPoints = zeros(1,numPoints);
            ListOfForbiddenPoints(1) = iidxPoint;

            % start of sorting: take the initial point and its nearest
            % neighbour
            IdxPointsSurfOrdered{itry}(2) = IidxNearestNeighbour(2, iidxPoint);
            ListOfForbiddenPoints(2) = IidxNearestNeighbour(2, iidxPoint);
            iidxPoint = IidxNearestNeighbour(2, iidxPoint);
            numPointsPerPart(itry) = 2;

            % consider the still unsorted points
            ipoint = 2;
            PointFound = true;

            while ipoint <= numPoints && PointFound
                
                % For the current point in the point set, consider the
                % nearest neighbour, if still available, otherwise the
                % second nearest neighbour, if still available and so on.
                % If this point passes the consistency check, we have found
                % the successor of the current point in ordering.
                % If we tried all knearestNeighbours next points and still
                % did not find a successor, we assume that all points of
                % the same segment are found and we proceed to the next
                % segment.
                PointFound = false;
                i = 1;
                while (i <= min(numPoints-1,knearestNeighbours) && ~PointFound)

                    % i-nearest neighbour is not forbidden
                    if (~any(ListOfForbiddenPoints == IidxNearestNeighbour(i+1, iidxPoint)))

                        % consistency check: compute the angle
                        % between the line segments: if it is smaller than
                        % alphaMax, we take the candidate and leave the
                        % loop. Otherwise, we reject the nearest available
                        % neighbour and try the next nearest one
                        seg1 = PointSet(iidxPoint,:) - PointSet(IidxNearestNeighbour(i+1, iidxPoint),:);
                        seg2 = PointSet(IdxPointsSurfOrdered{itry}(ipoint-1),:) - PointSet(iidxPoint,:); 

                        angle = (seg1*seg2')/(norm(seg1)*norm(seg2));
                        if (angle > FaultApproxParams.cosAlphaMax)
                            IdxPointsSurfOrdered{itry}(ipoint+1) = IidxNearestNeighbour(i+1, iidxPoint);
                            iidxPoint = IidxNearestNeighbour(i+1, iidxPoint);
                            ListOfForbiddenPoints(ipoint+1) = iidxPoint;
                            PointFound = true;
                            numPointsPerPart(itry) = numPointsPerPart(itry)+1;
                        end
                    end
                    i = i+1;
                end
                ipoint = ipoint +1;
            end

            % if there are no points to sort left: Stop sorting
            if ~any(IdxPointsSurfOrdered{itry} ==0)
                IdxPointsSurfOrdered{itry} = IdxAux(IdxPointsSurfOrdered{itry});
                ncomps = itry;
                PointFound = true;
                break

            % there are points left to sort: reduce the point set to these
            % leftover points and proceed
            else
                % shorten the index array to the correct length
                IdxPointsSurfOrdered{itry} = IdxPointsSurfOrdered{itry}(1:numPointsPerPart(itry));

                % indices with respect to the current point set
                IdxPointsSurfOrderedLocal = IdxPointsSurfOrdered{itry};

                % indices with respect to the original point set
                IdxPointsSurfOrdered{itry} = IdxAux(IdxPointsSurfOrdered{itry});

                % idx contains the indices of the points already sorted with
                % respect to the original point set

                idxLocal = setdiff(1:numPoints, IdxPointsSurfOrderedLocal);

                % indices of the points already sorted
                DistMat = DistMat(idxLocal, idxLocal);
                PointSet = PointSet(idxLocal,:);
                [~, IidxNearestNeighbour] = sort(DistMat);
                numPoints = size(PointSet,1);
                IdxAux = IdxAux(idxLocal);
            end
            
        % one orphaned point left: There is no need for sorting. Take the
        % one and only point and make it an own component
        else
            IdxPointsSurfOrdered{itry} = IdxAux(1);
            numPointsPerPart(itry) = 1;
            ncomps = itry;
            PointFound = true;
            break
        end
    end

    % There are still points which could not be included into the set of
    % sorted points and we have used all our trials: we failed.
    if (~PointFound && itry == ntrials)
        sortingSuccessful = false;
        
    % sorting was successful so far: try to combine the subsets if
    % desired and if there are more than one subsets
    else
        if bforceCombine && ncomps > 1

            % index of the nearest subset
            isubMin = 1;
            
            % distance between the starting points of the different
            % components
            normStartStart = zeros(1, ncomps);

            normStartEnd = zeros(1, ncomps);
            normEndStart = zeros(1, ncomps);

            % distance between the end points of the different
            % components
            normEndEnd = zeros(1, ncomps);

            % compute the distance between the starting and end points of
            % the different subsets
            for isub = 1 : ncomps-1
                
                % start and end points of the isub-th subset
                xstart = PointSetOrig(IdxPointsSurfOrdered{1}(1),:); 
                xend = PointSetOrig(IdxPointsSurfOrdered{1}(numPointsPerPart(1)),:);

                auxOld = 1e20;
                for jsub = 2: ncomps-isub+1

                    % find the segment which is closest.
                    startNext = PointSetOrig(IdxPointsSurfOrdered{jsub}(1),:);
                    endNext = PointSetOrig(IdxPointsSurfOrdered{jsub}(numPointsPerPart(jsub)),:);
                    normStartStart(jsub) = norm(xstart - startNext);
                    normStartEnd(jsub) = norm(xstart - endNext);
                    normEndStart(jsub) = norm(xend - startNext);
                    normEndEnd(jsub) = norm(endNext - xend);
                    aux = min([normStartStart(jsub), normStartEnd(jsub), normEndStart(jsub), normEndEnd(jsub)]);
                    if (aux < auxOld)
                        auxOld = aux;
                        isubMin = jsub;
                    end
                end

                % sew the two subsets together
                % the two starting points are closest
                mergeSuccess = false;
                if (normStartStart(isubMin) <= normStartEnd(isubMin) && ...
                    normStartStart(isubMin) <= normEndStart(isubMin) && ...
                    normStartStart(isubMin) <= normEndEnd(isubMin))
                    IdxPointsSurfOrdered{1} = [flip(IdxPointsSurfOrdered{1}) IdxPointsSurfOrdered{isubMin}];
                    mergeSuccess = true;
                % the starting point of the first and the end points of the
                % second subset are closest
                elseif (normStartEnd(isubMin) <= normStartStart(isubMin) && ...
                        normStartEnd(isubMin) <= normEndStart(isubMin) && ...
                        normStartEnd(isubMin) <= normEndEnd(isubMin))
                    IdxPointsSurfOrdered{1} = [flip(IdxPointsSurfOrdered{1}) flip(IdxPointsSurfOrdered{isubMin})];
                    mergeSuccess = true;
                elseif(normEndStart(isubMin) <= normStartStart(isubMin) && ...
                       normEndStart(isubMin) <= normStartEnd(isubMin) && ...
                       normEndStart(isubMin) <= normEndEnd(isubMin))
                    IdxPointsSurfOrdered{1} = [IdxPointsSurfOrdered{1} IdxPointsSurfOrdered{isubMin}];
                    mergeSuccess = true;
                elseif(normEndEnd(isubMin) <= normStartStart(isubMin) && ...
                       normEndEnd(isubMin) <= normStartEnd(isubMin) && ...
                       normEndEnd(isubMin) <= normEndStart(isubMin))
                    IdxPointsSurfOrdered{1} = [IdxPointsSurfOrdered{1} flip(IdxPointsSurfOrdered{isubMin})];
                    mergeSuccess = true;
                else
                    warning('Merging of parts failed in sortPointsOnFaultLine.')
                end

                % delete the subset just attached and continue
                % combine all subsets in the first one
                if mergeSuccess
                    numPointsPerPart(1) = numPointsPerPart(1) + numPointsPerPart(isubMin);
                    for i = isubMin:ncomps-1
                        IdxPointsSurfOrdered{i} = IdxPointsSurfOrdered{i+1};
                        numPointsPerPart(i) = numPointsPerPart(i+1);
                    end
                end                
            end
        end
        
        if (ncomps > 0)
            sortingSuccessful = true;
        end
    end
    
    % it may happen that the polygonal line after search intersects, but
    % only if the fault in fact consists of several components. Then,
    % setting bforceCombine = true sews components together in error
    % which do not belong together. However, if combining the subsets is
    % not forced, this must not happen.
    if (~bforceCombine)
        for icomp = 1:ncomps
            doesNotIntersect = ~selfIntersection(PointSetOrig(IdxPointsSurfOrdered{icomp},:));
            if (~doesNotIntersect)
                warning('Sorting leads to a self-intersecting polygonal line.')
            end
        end
    end    
    
    if (~sortingSuccessful)
        warning('Sorting of points on fault line failed.')
    end
end