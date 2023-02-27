% This function adaptively refines a set of triplets representing a fault
% surface based upon the estimated l_oo-error of a hypothetical
% triangulation using the given triplets.
% In three dimensions, there is no adaptive coarsening, as usually the
% resulting set of triplets is used for obtaining an approximative
% triangulation of the fault surface. However, the "holes" created by
% adaptive coarsening make meshing much more complicated.
%
% Input:
% - PointsIclass: Array containing the points which belong to class iclass.
%   Size: (number of points x 3)
% - PointsJclass: Array containing the points which belong to class jclass.
%   The arrays PointsIclass and PointsJclass are organized such that
%   PointsIclass(ipoint) and PointsJclass(ipoint) belong to the same
%   triplet.
%   Size: (number of points x 3)
% - iclass, jclass: class indices
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
%
% Output:
% - PointsIclass: Array containing the points which belong to class iclass.
%   Size: (number of points x 3)
% - PointsJclass: Array containing the points which belong to class jclass.
%   Size: (number of points x 3)
% - bsuccessful: boolean (true, if adaptive refinement and coarsening was
%   sucessful, false otherwise)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointsIclass, PointsJclass, bsuccessful] = adapt3D(PointsIclass, ...
    PointsJclass, iclass, jclass, FaultApproxParams, ProblemDescr)


    global ncalls

    bsuccessful = true;
    
    [npoints, ndim] = size(PointsIclass);
   
    % local tolerance
    epsLoc = 1e-6;
    
    % precomputed square for less cluttered code
    minDistFactorSq = FaultApproxParams.minDistFactor^2;
    
    % Flag which indicates if any points hve been added or not.
    refinementDone = false;
    
    for itry = 1:FaultApproxParams.maxiterAdapt
        % For a meaningful triangulation-based algorithm, we need at least
        % three points.
        if (npoints > 2)
            % Compute the distance matrix (actually the squares of
            % distances).
            DistMat = zeros(npoints, npoints, ndim);
            for k = 1: npoints
                DistMat(:,k,:) = PointsIclass-PointsIclass(k,:);
            end
            
            DistMat = DistMat.^2;
            for idim = 2: ndim
                DistMat(:,:,1) = DistMat(:,:,1) + DistMat(:,:,idim);
            end
            [DistMatSorted, IidxNearestNeighbour] = sort(DistMat(:,:,1),1);

            NewPoints = zeros(0,3);
            NewPointsNormals = zeros(0,3);
            NewStepSizes = zeros(0,1);
            NewDensity = zeros(0,1);

            avgDensity = zeros(npoints, 1);
            
            npointsAux = min(30, npoints);
            for ipoint = 1: npoints

                radiusSquared = DistMatSorted(npointsAux, ipoint);
                avgDensity(ipoint) = npointsAux/(pi*radiusSquared);
            end
            
            % When starting with a very coarse point set, npoints may even
            % smaller than nNearestPoints. However, after adding points,
            % this may not be the case any more. So we recompute
            % nNearestPoints for any filling iteration.
            nNearestPoints = min(FaultApproxParams.nNearestPoints, npoints);

            for ipoint = 1:npoints
                NearestPoints = 0.5*(PointsIclass(IidxNearestNeighbour(1:nNearestPoints, ipoint),:) + ...
                                     PointsJclass(IidxNearestNeighbour(1:nNearestPoints, ipoint),:));

                % midpoint aka center of gravity
                xmid = 1/nNearestPoints*ones(1, nNearestPoints)*NearestPoints;

                % compute local coordinate system
                % We compute the optimally fitting plane in the sense that
                % the sum of the squared distances between the points and
                % the plane is minimal.
                % It is well known (see e.g. Shakarji, M.: Least-Squares
                % Fitting Algorithms of the NIST Algorithm Testing System,
                % J. Res. Natl. Inst. Stand. Technol. 103, 633 (1998),
                % https://nvlpubs.nist.gov/nistpubs/jres/103/6/j36sha.pdf)
                % that the optimal fitting plane contains the midpoint.
                % Therefore, we shift the point such that 0 is the new
                % midpoint beforehand.
                NearestPointsShifted = NearestPoints - xmid;
                [Q, ~] = eig(NearestPointsShifted'*NearestPointsShifted);

                % "flatten" the neighbouring points by projecting them into
                % the optimally fitting plane.
                NearestPointsInPlane = NearestPointsShifted*Q;

                PointValues = NearestPointsInPlane(:,1);
                NearestPointsInPlane = NearestPointsInPlane(:,2:3);

                % Remove duplicate or almost duplicate points in the local
                % plane
                % Vectorisation makes this nested loop even slower!
                for i = 1: nNearestPoints - 1
                    for j = i+1: nNearestPoints
                        if(norm(NearestPointsInPlane(i,:) - NearestPointsInPlane(j,:)) < epsLoc)
                            NearestPointsInPlane(j,:) = -42;
                       end
                    end
                end
                
                NearestPoints = NearestPoints(NearestPointsInPlane(:,1) ~= -42,:);
                PointValues = PointValues(NearestPointsInPlane(:,1) ~= -42,:);
                NearestPointsInPlane = NearestPointsInPlane(NearestPointsInPlane(:,1) ~= -42,:);
                                
                if (size(NearestPointsInPlane,1) >= 3)
                    % Compute a 2d-triangulation using the Matlab-function
                    % delaunay.
                    ConnectivityList = delaunay(NearestPointsInPlane);

                    % Fill holes in triangulation.
                    % Project all nearest neighbours on the plane with
                    % normal vector NormalVec which contains point ipoint.
                    numberOfTriangles = size(ConnectivityList,1);
                    NewPointsLoc = zeros(numberOfTriangles,3);
                    NewPointsNormalsLoc = zeros(numberOfTriangles,3);
                    NewPointsTrisLoc = zeros(numberOfTriangles,3);
                    NewStepSizesLoc = zeros(numberOfTriangles,1);
                    NewDensityLoc = zeros(numberOfTriangles,1);
                    j = 1;
                    
                    edgeLengthAux = zeros(3,3);

                    [coeffs, scale] = compApproxRBF3D(NearestPointsInPlane, PointValues, 0.001);     
                    
                    for itri = 1: numberOfTriangles

                        % get local indices of the triangle vertices
                        iidx = ConnectivityList(itri,:);

                        if (any(iidx==1))
                        
                            edge1 = NearestPoints(iidx(1),:) - NearestPoints(iidx(2),:);
                            edge2 = NearestPoints(iidx(2),:) - NearestPoints(iidx(3),:);
                            edgeLengthAux(3,:) = NearestPoints(iidx(3),:) - NearestPoints(iidx(1),:);

                            % compute edge length on the real surface                    
                            edgeLengthAux(1,:) = edge1;
                            edgeLengthAux(2,:) = edge2;

                            % actually, the square of the maximal edge
                            % length
                            edgeLength = max((edgeLengthAux.^2)*[1;1;1]);

                            % Rough estimation of the normal vector. The
                            % function computeAllNormalVecs uses a way more
                            % sophisticated method. However, we need the
                            % normal just for computing a starting point
                            % for bisection. So, the cross product is just
                            % fine here.
                            % Moreover, we want to estimate the "flatness"
                            % of the triangle, which we do by estimating
                            % the ratio of area vs. edgeLength. The area of
                            % the triangle is || e1 x e2|| = norm(normal).
                            normal = cross(edge1, edge2);

                            % We refine the triangle if its edges are not
                            % too short and if it is not extremely flat:
                            % Consider these points in the plane with
                            % 2 = (1,1), 3 = (1,0.5), 5 = (1,0):
                            %  1--------------2
                            %  |              |   
                            %  |              |
                            %  |              3
                            %  |              |
                            %  |              |
                            %  4--------------5
                            % Then, the Delaunay triangulation is [(3,4,5),
                            % (1,3,2), (1,4,3)]. However, sometimes, the
                            % triangle (2,5,3) occurs, seemingly due to
                            % rounding errors. Therefore, we have to
                            % exclude this corner case.
                            if (norm(normal)/edgeLength > 0.01)

                                errMax = estimateMaxErr3D(NearestPointsInPlane, ...
                                    ConnectivityList, itri, coeffs, scale);                                  
                                
                                % We refine, if either the triangle is too
                                % large or the error is too large. As error
                                % estimation is reliable only if no aux
                                % points have been added.
                                if (errMax > FaultApproxParams.errMax)

                                    % global indices of the vertices
                                    % forming the current triangle
                                    iidxGlobal = [IidxNearestNeighbour(iidx(1), ipoint), ...
                                                  IidxNearestNeighbour(iidx(2), ipoint), ...
                                                  IidxNearestNeighbour(iidx(3), ipoint)];

                                    % As candidate for a new surface point,
                                    % take the barycentre of the triangle,
                                    % if appropriate.
                                    NewPointTest = 1/3*(NearestPointsInPlane(iidx(1),:) + ...
                                                        NearestPointsInPlane(iidx(2),:) + ...
                                                        NearestPointsInPlane(iidx(3),:));

                                    phi = 0;
                                    numPoints = size(NearestPointsInPlane, 1);
                                    ScaleVec = scale*ones(1, numPoints);
                                    for i = 1: numPoints
                                       phi = phi + coeffs(i)* Gaussian(NewPointTest, NearestPointsInPlane(i,:), ...
                                              ScaleVec(i));
                                    end

                                    NewPointsLoc(j,:) = [phi, NewPointTest]*Q' + xmid;

                                    NewPointsNormalsLoc(j,:) = normal;
                                    NewStepSizesLoc(j) = errMax;
                                    NewPointsTrisLoc(j,:) = iidxGlobal;
                                    NewDensityLoc(j) = 1/3*(avgDensity(iidx(1)) + ...
                                                            avgDensity(iidx(2)) + ...
                                                            avgDensity(iidx(3)));
                                    j = j + 1;
                                end
                            end
                        end
                    end
                    NewPoints = [NewPoints; NewPointsLoc(1:j-1,:)];
                    NewPointsNormals = [NewPointsNormals; NewPointsNormalsLoc(1:j-1,:)];
                    NewStepSizes = [NewStepSizes; NewStepSizesLoc(1:j-1)];
                    NewDensity = [NewDensity; NewDensityLoc(1:j-1)];

                % Actually, if this else-branch is entered, something went
                % completely wrong beforehand.
                else
                    warning('Local triangulation failed because of too less points.')
                    warning('This indicates that there are duplicates in the point set.')
                end
            end
            
            npointsNewGlob = size(NewPoints,1);
            minDist = FaultApproxParams.minDistFactor*FaultApproxParams.maxDistForSurfacePoints;

            if (npointsNewGlob == 0)
                refinementDone = true;
                break
            end
            
            % Test, if some of the new candidates are very close together
            % or even identical. This may happen as every points belongs to
            % several patches. It is possible that by some coincidence,
            % several patches consist the very same triangle. If this is
            % marked for refinement, the same points is added several
            % times.
            for i = 1: npointsNewGlob

                % square of distance of the i-the point to all subsequent
                % ones
                distVec = NewPoints(i,:) - NewPoints(i+1:end,:);
                distVec = (distVec.^2)*[1; 1; 1];
                
                % mark the subsequent points that are too close to the
                % i-the one
                jvec = false(npointsNewGlob,1);
                jvec(i+1:end) = distVec < min ( minDist*minDist, minDistFactorSq/NewDensity(i));
            
                % mark the almost duplicate points for removal
                NewPoints(jvec,:) = -42;
                NewPointsNormals(jvec,:) = -42;
                NewStepSizes(i) = NewStepSizes(i) + sum(NewStepSizes(jvec & NewStepSizes > -24))/ ...
                (sum(jvec & NewStepSizes > -24) + 1);
                NewStepSizes(jvec) = -42;
            end
                        
            % remove the almost duplicate points
            NewPoints = NewPoints(NewPoints(:,1) ~= -42,:);
            NewStepSizes = NewStepSizes(NewStepSizes ~= -42);
            NewPointsNormals = NewPointsNormals(NewPointsNormals(:,1) ~= -42,:);
            
            % Even if there are duplicates, at least one point must
            % remain. So it is pointless to break at this place.
            
            % Same test with respect to existing points. This is a very
            % heuristic approach: If a starting point is good and very
            % close to an existing point, it is likely that the final new
            % point is very close to an existing point as well. We want to
            % avoid building clusters. Therefore, we remove such
            % candidates.
            for i = 1: npoints
                distVec = 0.5*(PointsIclass(i,:) + PointsJclass(i,:)) - NewPoints(:,:);
                distVec = (distVec.^2)*[1; 1; 1];
                jvec = distVec < min(minDist*minDist, minDistFactorSq/avgDensity(i));
                NewPoints(jvec,:) = -42;
                NewPointsNormals(jvec,:) = -42;
                NewStepSizes(jvec) = -42;
            end
            NewPoints = NewPoints(NewPoints(:,1) ~= -42,:);
            NewStepSizes = NewStepSizes(NewStepSizes ~= -42);
            NewPointsNormals = NewPointsNormals(NewPointsNormals(:,1) ~= -42,:);

            % ensure all candidates are inside the domain
            NewStepSizes = NewStepSizes(all(NewPoints> ProblemDescr.Xmin,2) & ...
                                        all(NewPoints < ProblemDescr.Xmax,2));
            NewPointsNormals = NewPointsNormals(all(NewPoints> ProblemDescr.Xmin,2) &...
                                                all(NewPoints < ProblemDescr.Xmax,2),:);
            NewPoints = NewPoints(all(NewPoints> ProblemDescr.Xmin,2) & ...
                                  all(NewPoints < ProblemDescr.Xmax,2),:);

            npointsNewGlob = size(NewPoints,1);

            % If there are no new points left: all gaps must have been
            % filled.
            if (npointsNewGlob == 0)
                refinementDone = true;
                break
            end
            
            % normalise NewPointsNormals
            for i = 1: npointsNewGlob
                NewPointsNormals(i,:) = NewPointsNormals(i,:)/norm(NewPointsNormals(i,:));
            end
            
            %

            NewStepSizes = max(0.95*FaultApproxParams.abstolBisection, ...
                               0.5*min(0.25*NewStepSizes, ...
                                       FaultApproxParams.alpha*FaultApproxParams.maxDistForSurfacePoints));
            NewPointsBesidesSurface = NewPoints + NewPointsNormals.*kron(NewStepSizes, ones(1,3));
            NewPointsBesidesSurface = min(max(ProblemDescr.Xmin + eps, NewPointsBesidesSurface), ...
                                              ProblemDescr.Xmax - eps);

            NewPoints = NewPoints - NewPointsNormals.*kron(NewStepSizes, ones(1,3));
            NewPoints = min(max(ProblemDescr.Xmin + eps, NewPoints), ProblemDescr.Xmax - eps);


            AuxArr = computeClassification([NewPoints; NewPointsBesidesSurface], ProblemDescr);

            ClassPointsOnCurve = AuxArr(1: size(NewPoints,1), :);
            ClassPointsBesidesCurve = AuxArr(size(NewPoints,1)+1:size(AuxArr, 1));

            % Find points near the fault line, each with a counterpart
            % in the opposite class.
            [PointPairs, IdxSucceeded, ClassPointsSucceeded] = startPairs(NewPoints, ...
                NewPointsBesidesSurface, ClassPointsOnCurve, ClassPointsBesidesCurve, ...
                iclass, jclass, ProblemDescr);

            npointsOld = npoints;
            
            % Compute points near the fault surface by bisection.
            for i = 1: size(IdxSucceeded,1)
                
                if (IdxSucceeded(i))
                    [PointLeft, PointRight, finished] = singleTripletByBisection(PointPairs(i,1:ndim), ...
                        PointPairs(i, ndim+1:2*ndim), ClassPointsSucceeded(i,1), ...
                        ClassPointsSucceeded(i,2), ProblemDescr, FaultApproxParams);
                    if finished
                        npoints = npoints + 1;

                        if (ClassPointsSucceeded(i,1) == iclass)
                            PointsIclass(npoints,:) = PointLeft;
                            PointsJclass(npoints,:) = PointRight;
                        else
                            PointsJclass(npoints,:) = PointLeft;
                            PointsIclass(npoints,:) = PointRight;
                        end
                    else
                        IdxSucceeded(i) = 0;
                        warning('computation of new point failed')
                    end
                end
            end

            if (npoints == npointsOld)
                warning('Premature end of adaptive refinement, as all candidates proved unsuitable.')
                break
            end
                
            % NewPoints are on a triangle and not on the actual surface.
            % Therefore, it may happen, that a current point, not very
            % close to an existing point leads after bisection to an
            % existing point, which would be duplicate then. This is
            % because the new point has been a starting point for bisection
            % in an earlier iteration. However, we do
            % not keep in mind the NewPoint-vector from the previous
            % iteration, such that we cannot detect this corner case here.
            % This is why we need to classify first and to remove
            % turned-out duplicates afterwards.
            auxVec = PointsIclass(1:npointsOld,:);
            for i = npointsOld+1:npoints
                
                % distance of the i-th new point to all old points
                distVec = PointsIclass(i,:) - auxVec;
                distVec = (distVec.^2)*[1; 1; 1];

                jvec = distVec < epsLoc*epsLoc;
                auxVec(jvec,:) = -42;
                NewPointsNormals(jvec,:) = -42;
            end
            PointsIclass = [auxVec(auxVec(:,1) ~= -42,:); PointsIclass(npointsOld+1:npoints,:)];
            PointsJclass = [PointsJclass(auxVec(:,1) ~= -42,:); PointsJclass(npointsOld+1:npoints,:)];

            npoints = size(PointsIclass, 1);
            
            if (npoints == npointsOld)
                warning('All added points already existed in the point set.')
                break
            end
            
        % just two points or less
        else
            warning('Too less points available for adaptive refinement.')
        end
    end
    
    if (~refinementDone)
        warning('Adaptive refinement of the fault surface terminated, but not finished.')
    end
end