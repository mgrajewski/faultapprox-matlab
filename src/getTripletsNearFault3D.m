% This function provides the 3d-specific part of the algorithm. For
% details, we refer to "Detecting and approximating decision boundaries in
% low dimensional spaces".
%
% Input:
% - PointSetsSurface: (nclasses x nclasses)-structure of structures
%   containing the points describing the faults
% - NumPointsSurf (nclasses x nclasses)-structure of arrays containing
%   the number of points in the point sets
% - nclasses: total number of classes
% - ClassVals: Array containing the class values. These are not
%   necessarily the class indices. Imagine that f(\Omega) = {1,2,5}. Then,
%   ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
%   Size: nclasses
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
%
% Output:
% - PointSetsSurface: (nclasses x nclasses)-structure of structures
%   containing the points describing the faults
% - NormalsSurface: the same as PointSetsSurface, but with outer normal
%   vectors
% - NumPointsSurf (nclasses x nclasses)-structure of arrays containing
%   the number of points in the point sets
% - bsuccessful: flag, if the fault lines have been processed successfully

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointSetsSurface, NormalsSurface, NumPointsSurf, bsuccessful] = ...
    getTripletsNearFault3D(PointSetsSurface, NumPointsSurf, nclasses, ...
                           ClassVals, FaultApproxParams, ProblemDescr)
    NormalsSurface = cell(nclasses);

    global ncalls
    global ExtendedStats;

    alphaMaxExpand = 0.4;
    alphaExpandMin = 0.7;

    % points closer as epsLoc are regarded as identical
    epsLoc = FaultApproxParams.eps;

    % find out dimension
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses
            PointSet = PointSetsSurface{iclass, jclass};
            if (~isempty(PointSet))
                ndim = size(PointSet{1},2);
                break;
            end
        end
    end


    % we preallocate the arrays in order to avoid dynamic resizing. We
    % shorten them later.
    for iclass = 1:nclasses
        for jclass = 1:iclass -1
            NormalsSurface{iclass,jclass} = cell(1,1);
        end

        for jclass = iclass+1:nclasses
            NormalsSurface{iclass,jclass} = cell(1,1);
        end
    end

    bsuccessful = false;
    
    % indices of the domain boundaries
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
    % The expansion step is more complicated than in 2D: When the gaps
    % are filled, we compute the bounding box around the cloud of
    % points near the decision surface. 
    LeftDomainStart = zeros(nclasses, nclasses, 6);        
    LeftDomainEnd = zeros(nclasses, nclasses, 6);        
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses

            if ProblemDescr.verboseMode
                disp(['-- fill gaps on boundary between classes ', ...
                      int2str(ClassVals(iclass)), ' and ', ...
                      int2str(ClassVals(jclass))])
            end

            if (NumPointsSurf{iclass, jclass}(1) > 0)
                 [PointSetsSurface{iclass, jclass}{1}, PointSetsSurface{jclass, iclass}{1}, ...
                     successFill] = ...
                     fill3D(PointSetsSurface{iclass, jclass}{1}, PointSetsSurface{jclass, iclass}{1}, ...
                                            iclass, jclass, FaultApproxParams, ProblemDescr);

                NumPointsSurf{iclass, jclass}(1) = size(PointSetsSurface{iclass, jclass}{1},1);                    
                NumPointsSurf{jclass, iclass}(1) = NumPointsSurf{iclass, jclass}(1);
            end
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_filling_gaps';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end


    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses

            if ProblemDescr.verboseMode
                disp(['-- expand boundary between classes ' ...
                      int2str(ClassVals(iclass)) ' and ' ...
                      int2str(ClassVals(jclass))])
            end

            if (NumPointsSurf{iclass, jclass}(1) > 0)

                % expand the surface until its boundaries
                
                % compute the size of the bounding box
                BoundingBox = zeros(3,1);
                BoundingBox(1) = max(PointSetsSurface{iclass, jclass}{1}(:,1)) - ...
                                 min(PointSetsSurface{iclass, jclass}{1}(:,1));
                BoundingBox(2) = max(PointSetsSurface{iclass, jclass}{1}(:,2)) - ...
                                 min(PointSetsSurface{iclass, jclass}{1}(:,2));
                BoundingBox(3) = max(PointSetsSurface{iclass, jclass}{1}(:,3)) - ...
                                 min(PointSetsSurface{iclass, jclass}{1}(:,3));

                PointsNearBdryIclass = zeros(0);
                PointsNearBdryJclass = zeros(0);
                OuterNormalsDomain = zeros(0);
                iidxMax = zeros(0);
                
                % store all new points on the boundary in separate arrays
                NewPointsOnBoundaryIclass = cell(1,6+nclasses);
                NewPointsOnBoundaryJclass = cell(1,6+nclasses);

                Normals3D = computeAllNormalVecs(PointSetsSurface{iclass, jclass}{1}, ...
                                                 PointSetsSurface{jclass, iclass}{1}, ...
                                                 FaultApproxParams);

                % expansion to inner intersection lines
                for iclassTest = 1:nclasses
                    for jclassTest = iclassTest+1:nclasses
                        if (iclassTest ~= iclass && NumPointsSurf{iclassTest, jclassTest}(1) > 0)

                            % compute the distance matrix between the
                            % points in the current PointSet and the
                                % another one
                            DistMatAux = zeros(NumPointsSurf{iclass, jclass}(1), ...
                                               NumPointsSurf{iclassTest, jclassTest}(1), ...
                                               ndim);
                            for k = 1: NumPointsSurf{iclassTest, jclassTest}(1)
                                DistMatAux(:,k,:) = PointSetsSurface{iclass,jclass}{1} - ...
                                                    PointSetsSurface{iclassTest,jclassTest}{1}(k,:);
                            end
            
                            DistMatAux = DistMatAux.^2;
                            for idim = 2: ndim
                                DistMatAux(:,:,1) = DistMatAux(:,:,1) + DistMatAux(:,:,idim);
                            end
                            DistMat = DistMatAux(:,:,1).^0.5;
                            
                            % find out, if there are points sufficiently
                            % close to the other fault surface aka
                            % closer than maxDistForSurfacePoints
                            [DistToNearestPoint, minIndex] = min(DistMat, [], 2);

                            CloseToInnerBdry = DistToNearestPoint < 0.75*FaultApproxParams.maxDistForSurfacePoints;
                            PointsCloseToInnerBdry = PointSetsSurface{iclass,jclass}{1}(CloseToInnerBdry, :);
                            PointsNearBdryIclass = [PointsNearBdryIclass; PointsCloseToInnerBdry];
                            PointsNearBdryJclass = [PointsNearBdryJclass; PointSetsSurface{jclass,iclass}{1}(CloseToInnerBdry, :)];
                            
                            % get the approximate normal vector of surface
                            NormalsToSurface = Normals3D(CloseToInnerBdry,:);

                            % As we want to expand the surface along a
                            % line, we need a suitable direction. For
                            % doing so, we compute the vector to the
                            % nearest point on the other fault surface
                            % and project it into the tangential plane
                            % given by the estimated normal vector
                            NormalsInSurface = ...
                                PointsCloseToInnerBdry - PointSetsSurface{iclassTest,jclassTest}{1}(minIndex(CloseToInnerBdry),:);

                            % However, we need a vector perpendicular to
                            % the direction of the line to expand along
                            % and tangential to the surface, so we
                            % compute the cross product with
                            % NormalsToSurface
                            for i = 1: size(NormalsInSurface, 1)
                                NormalsInSurface(i,:) = NormalsInSurface(i,:) - ...
                                    [NormalsInSurface(i,:)*NormalsToSurface(i,:)'*NormalsToSurface(i,:)];
                                NormalsInSurface(i,:) = cross(NormalsInSurface(i,:), NormalsToSurface(i,:));
                            end
                            NormalsInSurface = -NormalsInSurface./vecnorm(NormalsInSurface, 2,2);
                            OuterNormalsDomain = [OuterNormalsDomain; NormalsInSurface];
                            iidxMax = [iidxMax; find(CloseToInnerBdry ~= 0)];
                       end

                    end
                end

                % expansion to the domain boundary
                for idim = 1: 3
                    
                    % identify candidates for expansion based upon
                    % the relative size of the elongation with
                    % respect to the idim-th and idim-1st coordinate
                    npoints = max(BoundingBox(mod(idim,3)+1), BoundingBox(mod(idim+1,3)+1))/ ...
                        FaultApproxParams.maxDistForSurfacePoints;
                    npoints = ceil(npoints);

                    % find points with maximal idim-th coordinates
                    [~,iidxAux] = maxk(PointSetsSurface{iclass, jclass}{1}(:,idim),npoints);

                    PointsNearBdryIclass = [PointsNearBdryIclass;PointSetsSurface{iclass, jclass}{1}(iidxAux,:)]; %#ok<AGROW>
                    PointsNearBdryJclass = [PointsNearBdryJclass;PointSetsSurface{jclass, iclass}{1}(iidxAux,:)]; %#ok<AGROW>

                    OuterNormalsDomain = [OuterNormalsDomain; ones(npoints,1)*[idim ==1, idim ==2, idim==3]]; %#ok<AGROW>
                    iidxMax = [iidxMax; iidxAux]; %#ok<AGROW>
                    
                    % find points with minimal idim-th coordinates
                    [~,iidxAux] = mink(PointSetsSurface{iclass, jclass}{1}(:,idim),npoints);
                    PointsNearBdryIclass = [PointsNearBdryIclass; PointSetsSurface{iclass, jclass}{1}(iidxAux,:)]; %#ok<AGROW>
                    PointsNearBdryJclass = [PointsNearBdryJclass; PointSetsSurface{jclass, iclass}{1}(iidxAux,:)]; %#ok<AGROW>
                    OuterNormalsDomain = [OuterNormalsDomain; -ones(npoints,1)*[idim ==1, idim ==2, idim==3]]; %#ok<AGROW>
                    iidxMax = [iidxMax; iidxAux]; %#ok<AGROW>
                end

                npoints = size(PointsNearBdryIclass,1);
    
                % remove duplicate points
                % inner loop is pointless if the point with index i has been
                % removed before
                for i = 1: npoints
                    if (PointsNearBdryIclass(i,1) ~= -42)
                        for j = i+1: npoints
                           if(abs(norm(PointsNearBdryIclass(i,:) - PointsNearBdryIclass(j,:))) < epsLoc)
                                PointsNearBdryIclass(j,:) = -42;
                                PointsNearBdryJclass(j,:) = -42;
                           end
                        end
                    end
                end
                OuterNormalsDomain = OuterNormalsDomain(PointsNearBdryIclass(:,1) ~= -42,:);
                iidxMax = iidxMax(PointsNearBdryIclass(:,1) ~= -42,:);
                PointsNearBdryIclass = PointsNearBdryIclass(PointsNearBdryIclass(:,1) ~= -42,:);
                PointsNearBdryJclass = PointsNearBdryJclass(PointsNearBdryJclass(:,1) ~= -42,:);

                npoints = size(PointsNearBdryIclass,1);

                % compute the distance matrix with respect to the selected points
                DistMatAux = zeros(NumPointsSurf{iclass, jclass}(1), npoints, ndim);
                for k = 1: npoints
                    DistMatAux(:,k,:) = PointSetsSurface{iclass,jclass}{1} - PointsNearBdryIclass(k,:);
                end

                DistMatAux = DistMatAux.^2;
                for idim = 2: ndim
                    DistMatAux(:,:,1) = DistMatAux(:,:,1) + DistMatAux(:,:,idim);
                end
                DistMat = DistMatAux(:,:,1).^0.5;

                [~, IidxNearestNeighbour] = sort(DistMat,1);

                NormalsToSurface = computeAllNormalVecs(PointSetsSurface{iclass, jclass}{1}, ...
                                                        PointSetsSurface{jclass, iclass}{1}, ...
                                                        FaultApproxParams);
                
                
                for ipoint = 1: npoints
                    % find or compute a point such that we can expand
                    % these two points along the boundary of the fault
                    % surface
                    NormalVec = NormalsToSurface(iidxMax(ipoint),:);
                    
                    % the following might happen (in 2D)
                    %
                    %    \ decision surface
                    %     x
                    %      \
                    %       |
                    %       x <- rightmost point, but not suitable for
                    %      /     expanding the surface in x-direction
                    %     x      (indicated by a normal vector pointing 
                    %    /        almost in x-direction)
                    %  _/
                    % /
                    %
                    % We exclude this case by considering the angle
                    % between the (estimated) normal vector and the
                    % normal vector of the coordinate plane aka
                    % OuterNormals.

                    % we do need the absolute value as the normal vector
                    % may point into the wrong direction anyway
                    angle = abs(OuterNormalsDomain(ipoint,:)*NormalVec');
                    
                    % If the angle is less than approx.
                    % acos(alphaMaxExpand), the case displayed above
                    % is not given here. However, if the points are really
                    % close to the boundary, it is likely that the surface
                    % extends to the domain boundary.
                    if angle < alphaMaxExpand || (any(PointsNearBdryIclass(ipoint,:) > ...
                                                  ProblemDescr.Xmax - FaultApproxParams.maxDistForSurfacePoints)) ...
                                              || (any(PointsNearBdryIclass(ipoint,:) < ...
                                                  ProblemDescr.Xmin + FaultApproxParams.maxDistForSurfacePoints))
                    
                        % find vector perpendicular to NormalVec in the plane
                        % spanned by NormalVec and OuterNormals by using
                        % orthogonal projection
                        % This is approximately the direction to ideally
                        % proceed with in expanding to the domain boundary.
                        NormalVecInSurface = OuterNormalsDomain(ipoint,:);
                        VecForExpansion = cross(NormalVec, NormalVecInSurface);

                        % for expanding the surface we proceed as follows:
                        % We start with two initial points and then we
                        % proceed like in the 2D-case aka in expand2D
                        PointsAuxIclass = zeros(2,3);
                        PointsAuxJclass = zeros(2,3);

                        PointsAuxIclass(2,:) = PointsNearBdryIclass(ipoint,:);
                        PointsAuxJclass(2,:) = PointsNearBdryJclass(ipoint,:);
                        bfound = false;
                        nNearestPoints = min(FaultApproxParams.nNearestPoints+1, npoints) -1;
                        VecToNearestPoints = PointSetsSurface{iclass,jclass}{1} ...
                            (IidxNearestNeighbour(2:nNearestPoints+1, ipoint),:) - PointsNearBdryIclass(ipoint,:);

                        % for expanding on the fault surface to the
                        % domain boundary, we need apart of the
                        % existing point another one which is rather
                        % close and points approximately in direction
                        % of VecToProceed. If we find such a point, we
                        % can try to expand.
                        for i = 1:nNearestPoints

                            % angle between -VecToNearestPoints and
                            % VecForExpansion is at most 25°
                            angleLoc = -VecToNearestPoints(i,:)*VecForExpansion' /(norm(VecForExpansion)*norm(VecToNearestPoints(i,:)));
                            if (angleLoc > alphaExpandMin && DistMat(IidxNearestNeighbour(i+1, ipoint),ipoint) < FaultApproxParams.maxDistForSurfacePoints)
                                PointsAuxIclass(1,:) = PointSetsSurface{iclass,jclass}{1}(IidxNearestNeighbour(i+1, ipoint),:);
                                PointsAuxJclass(1,:) = PointSetsSurface{jclass,iclass}{1}(IidxNearestNeighbour(i+1, ipoint),:);
                                bfound = true;
                                break
                            end
                        end

                        if (bfound)

                            IdxPointsSurfOrdered = [1 2];
                            LeftDomain = zeros(nclasses,nclasses,1);


                            % expand inside the fault surface to its
                            % boundary to get points on the boundary of the
                            % fault surface
                            [IdxPointsSurfOrderedEnlarged, LeftDomain, numPoints, PointsAuxIclass, PointsAuxJclass, resort, pointsAdded] = ...
                               expand2D(LeftDomain, IdxPointsSurfOrdered, iclass, jclass, 1, ClassVals, ...
                                                   PointsAuxIclass, PointsAuxJclass, VecForExpansion, ...
                                                   ProblemDescr, FaultApproxParams, 1);

                            iplane = LeftDomain(iclass,jclass);                                

                            % expansion to domain boundaries
                            if (iplane ~= 0 && pointsAdded)

                                succeeded = false;
                                % test, if the points found are really on
                                % the domain boundary if advertised.
                                switch iplane
                                case 1
                                    if PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),2) < ProblemDescr.Xmin(2) + 3*FaultApproxParams.abstolBisection
                                        PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),2) = ProblemDescr.Xmin(2) + eps;
                                        PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),2) = ProblemDescr.Xmin(2) + eps;
                                        succeeded = true;
                                    end
                                case 2
                                    if PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),1) > ProblemDescr.Xmax(1) - 3*FaultApproxParams.abstolBisection
                                        PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),1) = ProblemDescr.Xmax(1) - epsLoc;
                                        PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),1) = ProblemDescr.Xmax(1) - epsLoc;
                                        succeeded = true;
                                    end
                                case 3
                                    if PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),2) > ProblemDescr.Xmax(2) - 3*FaultApproxParams.abstolBisection
                                        PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),2) = ProblemDescr.Xmax(2) - epsLoc;
                                        PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),2) = ProblemDescr.Xmax(2) - epsLoc;
                                        succeeded = true;
                                    end
                                case 4
                                    if PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),1) < ProblemDescr.Xmin(1) + 3*FaultApproxParams.abstolBisection
                                        PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),1) = ProblemDescr.Xmin(1) + epsLoc;
                                        PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),1) = ProblemDescr.Xmin(1) + epsLoc;
                                        succeeded = true;
                                    end
                                case 5
                                    if PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),3) < ProblemDescr.Xmin(3) + 3*FaultApproxParams.abstolBisection
                                        PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),3) = ProblemDescr.Xmin(3) + epsLoc;
                                        PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),3) = ProblemDescr.Xmin(3) + epsLoc;
                                        succeeded = true;
                                    end
                                case 6
                                    if PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),3) > ProblemDescr.Xmax(3) - 3*FaultApproxParams.abstolBisection
                                        PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),3) = ProblemDescr.Xmax(3) - epsLoc;
                                        PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),3) = ProblemDescr.Xmax(3) - epsLoc;
                                        succeeded = true;
                                    end
                                end

                                if (succeeded)
                                    numBndPoints = size(NewPointsOnBoundaryIclass{iplane},1);
                                    NewPointsOnBoundaryIclass{iplane}(numBndPoints+1,:) = PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),:);
                                    NewPointsOnBoundaryJclass{iplane}(numBndPoints+1,:) = PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),:);
    
                                    PointSetsSurface{iclass, jclass}{1} = [PointSetsSurface{iclass, jclass}{1}; PointsAuxIclass(1:numPoints-1,:)];
                                    PointSetsSurface{jclass, iclass}{1} = [PointSetsSurface{jclass, iclass}{1}; PointsAuxJclass(1:numPoints-1,:)];
    
                                    NumPointsSurf{iclass, jclass} = NumPointsSurf{iclass,jclass} + numPoints-1;
                                    NumPointsSurf{jclass, iclass} = NumPointsSurf{jclass,iclass} + numPoints-1;
                                end

                            % expansion to inner intersection curve
                            elseif (iplane == 0 && pointsAdded)
                                numBndPoints = size(NewPointsOnBoundaryIclass{6+jclass},1);

                                NewPointsOnBoundaryIclass{6+jclass}(numBndPoints+1,:) = PointsAuxIclass(IdxPointsSurfOrderedEnlarged(end),:);
                                NewPointsOnBoundaryJclass{6+jclass}(numBndPoints+1,:) = PointsAuxJclass(IdxPointsSurfOrderedEnlarged(end),:);

                                PointSetsSurface{iclass, jclass}{1} = [PointSetsSurface{iclass, jclass}{1}; PointsAuxIclass(1:numPoints,:)];
                                PointSetsSurface{jclass, iclass}{1} = [PointSetsSurface{jclass, iclass}{1}; PointsAuxJclass(1:numPoints,:)];

                                NumPointsSurf{iclass, jclass} = NumPointsSurf{iclass,jclass} + numPoints;
                                NumPointsSurf{jclass, iclass} = NumPointsSurf{jclass,iclass} + numPoints;
                            end
                        end
                    end
                end
                % if the fault surfaces intersect the domain
                % boundaries, expand these intersection lines as
                % far as possible
                NormalVecPlane = [0,-1,0; 1 0 0; 0 1 0; -1 0 0; 0 0 -1; 0 0 1];
                for iplane = 1: 6
                    numPoints = size(NewPointsOnBoundaryIclass{iplane},1);
                    if numPoints > 1

                        % The points on the boundary of the
                        % fault surface form a curve consisting of
                        % several parts, as the points are on the
                        % intersection of the current fault surface
                        % with either another one or with one of the
                        % facets of the outer boundary.
                        % We treat each "outer" part as a 2D fault line.
                        [IdxPointsSurfOrdered, sortingSuccessful] = ...
                            sortPointsOnFaultLine(NewPointsOnBoundaryIclass{iplane}, ...
                                                  1, ProblemDescr, ...
                                                  FaultApproxParams);
                        NewPointsOnBoundaryIclass{iplane} = ...
                            NewPointsOnBoundaryIclass{iplane}(IdxPointsSurfOrdered{1},:);
                        NewPointsOnBoundaryJclass{iplane} = ...
                            NewPointsOnBoundaryJclass{iplane}(IdxPointsSurfOrdered{1},:);

                        if (~sortingSuccessful)
                            warning('Sorting of points failed when computing intersections with domain boundaries.')
                            bsuccessful = false;
                            
                            return
                        end
                        
                        [NewPointsOnBoundaryIclass{iplane}, NewPointsOnBoundaryJclass{iplane}] = ...
                            removeClusters(NewPointsOnBoundaryIclass{iplane}, ...
                                           NewPointsOnBoundaryJclass{iplane}, ...
                                           numPoints, FaultApproxParams);


                        pointsAdded = true;
                        itry = 1;
                        while pointsAdded && itry <= FaultApproxParams.maxTrialsForFillingGaps
                            [NewPointsOnBoundaryIclass{iplane}, ...
                             NewPointsOnBoundaryJclass{iplane}, ~, ...
                             pointsAdded, FillingSuccessful] = ...
                                fill2D(NewPointsOnBoundaryIclass{iplane}, ...
                                       NewPointsOnBoundaryJclass{iplane}, ...
                                       iclass, jclass, ClassVals, ...
                                       FaultApproxParams, ProblemDescr, ...
                                       NormalVecPlane(iplane,:));
                                           
                            if (~FillingSuccessful)
                                break
                            end
                            % sort again
                            [IdxPointsSurfOrdered, sortingSuccessful] = ...
                                sortPointsOnFaultLine(NewPointsOnBoundaryIclass{iplane}, ...
                                                      1, ProblemDescr, ...
                                                      FaultApproxParams);
                            IdxPointsSurfOrdered = IdxPointsSurfOrdered{1};
                            itry = itry+1;
                            if (~sortingSuccessful)
                                warning('Sorting of points failed when computing intersections with domain boundaries.')
                                bsuccessful = false;
    
                                % dummy return values
                                LeftDomainStart = 0;
                                LeftDomainEnd = 0;
                                
                                return
                            end
                        end
                        
                        
                        [IdxPointsSurfOrderedEnlarged, LeftDomain, numPoints, NewPointsOnBoundaryIclass{iplane}, ...
                            NewPointsOnBoundaryJclass{iplane}, resort] = ...
                            expand2D(LeftDomain, IdxPointsSurfOrdered, iclass, jclass, 1, ClassVals, ...
                                     NewPointsOnBoundaryIclass{iplane}, ...
                                     NewPointsOnBoundaryJclass{iplane}, NormalVecPlane(iplane,:), ...
                                     ProblemDescr, FaultApproxParams, 0);
                        
                        if (resort)
                            [IdxPointsSurfOrderedEnlarged, sortingSuccessful] = ...
                                sortPointsOnFaultLine(NewPointsOnBoundaryIclass{iplane}, 1, ProblemDescr, FaultApproxParams);
                            IdxPointsSurfOrderedEnlarged = IdxPointsSurfOrderedEnlarged{1};
                        end

                        IdxPointsSurfOrdered = IdxPointsSurfOrderedEnlarged;
                        [IdxPointsSurfOrderedEnlarged, LeftDomain, numPoints, NewPointsOnBoundaryIclass{iplane}, ...
                            NewPointsOnBoundaryJclass{iplane}, resort] = ...
                            expand2D(LeftDomain, IdxPointsSurfOrdered, iclass, jclass, 1, ClassVals, ...
                                     NewPointsOnBoundaryIclass{iplane}, ...
                                     NewPointsOnBoundaryJclass{iplane}, NormalVecPlane(iplane,:), ...
                                     ProblemDescr, FaultApproxParams, 1);

                        if (resort)
                            [IdxPointsSurfOrderedEnlarged, sortingSuccessful] = ...
                                sortPointsOnFaultLine(NewPointsOnBoundaryIclass{iplane}, 1, ProblemDescr, FaultApproxParams);
                            IdxPointsSurfOrderedEnlarged = IdxPointsSurfOrderedEnlarged{1};
                        end

                        NewPointsOnBoundaryIclass{iplane} = ...
                            NewPointsOnBoundaryIclass{iplane}(IdxPointsSurfOrderedEnlarged,:);
                        NewPointsOnBoundaryJclass{iplane} = ...
                            NewPointsOnBoundaryJclass{iplane}(IdxPointsSurfOrderedEnlarged,:);
                        
                        PointSetsSurface{iclass,jclass}{1} = [PointSetsSurface{iclass,jclass}{1}; ...
                                                              NewPointsOnBoundaryIclass{iplane}];
                        PointSetsSurface{jclass,iclass}{1} = [PointSetsSurface{jclass,iclass}{1}; ...
                                                              NewPointsOnBoundaryJclass{iplane}];
                        NumPointsSurf{iclass, jclass} = NumPointsSurf{iclass,jclass} + numPoints;
                        NumPointsSurf{jclass, iclass} = NumPointsSurf{jclass,iclass} + numPoints;
                    end
                end
                
                if ProblemDescr.verboseMode
                    disp(['-- fill gaps on extended boundary between classes ' ...
                          int2str(ClassVals(iclass)) ' and ' ...
                          int2str(ClassVals(jclass))])
                end
                [PointSetsSurface{iclass, jclass}{1}, ...
                 PointSetsSurface{jclass, iclass}{1}, bsuccessful] = ...
                     fill3D(PointSetsSurface{iclass, jclass}{1}, ...
                            PointSetsSurface{jclass, iclass}{1}, ...
                            iclass, jclass, FaultApproxParams, ...
                            ProblemDescr);
            end % end fault surface is not empty
        end
    end

    % during expansion, duplicate points are possibly produced. We
    % remove them afterwards here. 
    for iclass = 1: nclasses -1

        for jclass = iclass + 1: nclasses
            numAux = NumPointsSurf{iclass, jclass}(1);

            % inner loop is pointless if the point with index i has been
            % removed before
            for i = 1: numAux
                if (PointSetsSurface{iclass, jclass}{1}(i,1) ~= -42)
                    for j = i+1: numAux
                       if(abs(norm(PointSetsSurface{iclass, jclass}{1}(i,:) - ...
                                   PointSetsSurface{iclass, jclass}{1}(j,:))) < epsLoc)
                            PointSetsSurface{iclass, jclass}{1}(j,:) = -42;
                            PointSetsSurface{jclass, iclass}{1}(j,:) = -42;
                       end
                    end
                end
            end
            PointSetsSurface{iclass, jclass}{1} = ...
                PointSetsSurface{iclass, jclass}{1}(PointSetsSurface{iclass, jclass}{1}(:,1) ~= -42,:);
            PointSetsSurface{jclass, iclass}{1} = ...
                PointSetsSurface{jclass, iclass}{1}(PointSetsSurface{jclass, iclass}{1}(:,1) ~= -42,:);

            NumPointsSurf{iclass, jclass} = size(PointSetsSurface{iclass, jclass}{1}, 1);
            NumPointsSurf{jclass, iclass} = size(PointSetsSurface{jclass, iclass}{1}, 1);
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_expand_surf';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end

    % adaptive refinement
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses

            if (NumPointsSurf{iclass, jclass}(1) > 0)

                if ProblemDescr.verboseMode
                    disp(['-- adaptive refinement on boundary between classes ' ...
                          int2str(ClassVals(iclass)) ' and ' ...
                          int2str(ClassVals(jclass))])
                end
                [PointSetsSurface{iclass, jclass}{1}, PointSetsSurface{jclass, iclass}{1}, bsuccessful] = ...
                     adapt3D(PointSetsSurface{iclass, jclass}{1}, PointSetsSurface{jclass, iclass}{1}, ...
                     iclass, jclass, FaultApproxParams, ProblemDescr);

                 
                 NumPointsSurf{iclass, jclass}(1) = size(PointSetsSurface{iclass, jclass}{1},1);                    
                 NumPointsSurf{jclass, iclass}(1) = NumPointsSurf{iclass, jclass}(1);

                NormalsToSurface = computeAllNormalVecs(PointSetsSurface{iclass, jclass}{1}, ...
                                                        PointSetsSurface{jclass, iclass}{1}, ...
                                                        FaultApproxParams);
                NormalsSurface{iclass, jclass}{1} = NormalsToSurface;
            end 
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_adaptive_ref';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end
    bsuccessful = true;
end