% This function constructs points near the fault line(s)/surfaces.
% For details, we refer to "Detecting and approximating decision boundaries
% in low dimensional spaces". We sketch our algorithm in what follows:
% 1) building block initialize: We subdivide the set of barycentres
%    according to class and search for any barycentre in class iclass
%    its nearest neighbour in class jclass. The line from iclass to
%    jclass is supposed to intersect the fault line/surface separating
%    subset iclass from subset jclass. We find two points, one in class
%    iclass and the other one in class jclass on this line close to the
%    fault line by bisection.
% 2) building block fill:
%    2D: We order these points according to their position on the fault
%    line. This enables us to compute the distance to the next point on
%    the line in the same class. If the distance to that point is too
%    large, we try to add more points between these two. If after
%    adding such points, there are still significant gaps between
%    consecutive points, this indicates that the fault line between
%    class iclass and class jclass consists of several separated
%    components. In this case, we split the set of points on the
%    current fault line into groups where each group represents one
%    component.
%    3D: see explanation in 
% 3) building block expand:
%    We try to expand each component of the fault line after the last
%    and before the first known point on that component in 2D or try to
%    expand the decision surface until its boundaries in 3D.
% 4) building block adapt: We improve the approximation of the
%    components of the fault lines/surfaces by adding or removing
%    points according to (estimated) curvature.
% 
% Input:
% - barycentres: (m x d)-array containing the cartesian coordinates of
%   the barycentres; d: dimension, m: number of points
% - PointSet: (m x d)-array containing the cartesian coordinates of
%   the sampling points; d: dimension, m: number of sampling points
% - ClassOfBarys: class indices of the barycentres
% - ClassOfPoints: class indices of the sampling points in PointSet
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - PointSetsSurface: (nclasses x nclasses)-array of structures
%   PointSetsSurface{i,j} refers to point sets between
%   classes i and j which belong to class i. It is an array of
%   structures and not just arrays as the intersection of classes i
%   and j may consist of several components, which can be accessed by
%   PointSetsSurface{i,j}{icomp}.
% - LeftDomainStart: in 2D, this
%   (nclasses x nclasses x #components)-array codes information on
%   whether a fault line starts on the domain boundary (value > 0, we
%   refer to getTripletsNearFault2D.m for ASCII-art explaining how) or
%   inside the domain (value = 0).
% - LeftDomainEnd: The same as LeftDomainStart, but for ending of fault
%   lines
% - bsuccessful: flag, if the fault lines have been found successfully
% - NormalsSurface: in 2D, just a dummy return value, in 3D: it is a
%   (nclasses x nclasses)-array of structures containing approximate
%   outer normals for all points in several components of the point sets
%   (organized analogously to PointSetsSurface).

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointSetsSurface, LeftDomainStart, LeftDomainEnd, bsuccessful, NormalsSurface] = ...
    getTripletsNearFaults(Barycentres, PointSet, ClassOfBarys, ...
                          ClassOfPoints, ProblemDescr, FaultApproxParams)

    % initialise with dummy return values
    LeftDomainStart = 0;
    LeftDomainEnd = 0;

    bsuccessful = false;
    global ncalls
    global ExtendedStats;
    
    % number of barycentres
    nbarys = size(Barycentres, 1);

    % dimension
    ndim = size(PointSet, 2);
    
    if (ndim < 2 || ndim > 3)
        error('The dimension of the domain must be 2 or 3.')
    end

    % Array containing the class values. These are not necessarily the
    % class indices. Imagine that f(\Omega) = {1,2,5}. Then,
    % ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
    ClassVals = unique(ClassOfPoints);
    
    % number of different classes in the point set
    nclasses = size(ClassVals,1);
    
    if (any(ClassVals < 1))
        error('Class values < 1 occurred. Classes must be however natural numbers.')
    end
    
    NormalsSurface = 0;
    
    % Get the opposite mapping: Class value -> i
    ClassValsInv = zeros(1, max(ClassVals));
    ClassValsInv(ClassVals) = 1:nclasses;
        
    % Possibly, there will be several lines or surfaces separating the
    % different classes. We store the points near the surface as follows:
    % All points in class i that belong to the surface separating class i
    % from j are stored in PointSetsSurface{i,j}, all points in class j
    % that belong to the surface separating class i
    % from j are stored in PointSetsSurface{j,i}.
    PointSetsSurface = cell(nclasses);

    % number of points on fault line PointSetsSurface{i,j}
    NumPointsSurf = cell(nclasses, nclasses);

    % We preallocate the arrays in order to avoid dynamic resizing. We
    % shorten them later.
    for iclass = 1:nclasses
        for jclass = 1:iclass -1
            PointSetsSurface{iclass,jclass} = cell(1,1);
            PointSetsSurface{iclass,jclass}{1} = -42*ones(nbarys, ndim);
            NumPointsSurf{iclass, jclass} = 0;
        end

        NumPointsSurf{iclass, iclass} = 0;

        for jclass = iclass+1:nclasses
            PointSetsSurface{iclass,jclass} = cell(1,1);
            PointSetsSurface{iclass,jclass}{1} = -42*ones(nbarys, ndim);
            NumPointsSurf{iclass, jclass} = 0;
        end
    end
    
    for iclass = 1: nclasses

        idxIclass = ClassVals(iclass);

        if ProblemDescr.verboseMode
            disp(['-- compute initial set of boundary points, class ', ...
                  int2str(idxIclass)])
        end
        
        % all points and barycentres not in the current class
        PointsNotInClass = [PointSet(ClassOfPoints~= idxIclass,:); ...
                            Barycentres(ClassOfBarys ~= idxIclass, :)];

        ClassOfPointsAux = [ClassOfPoints(ClassOfPoints~= idxIclass); ...
                            ClassOfBarys(ClassOfBarys ~= idxIclass)];

        % all barycentres in current class
        BarycentresInClass = Barycentres(ClassOfBarys == idxIclass,:);
        
        nbarysInClass = size(BarycentresInClass, 1);
        npointsNotInClass =size(PointsNotInClass, 1);

        DistMatAux = zeros(npointsNotInClass, nbarysInClass, ndim);
        for jclass = 1: nbarysInClass
            DistMatAux(:,jclass,:) = PointsNotInClass-BarycentresInClass(jclass,:);
        end
    
        DistMatAux = DistMatAux.^2;
        for idim = 2: ndim
            DistMatAux(:,:,1) = DistMatAux(:,:,1) + DistMatAux(:,:,idim);
        end
        
        % distance matrix to find the nearest point to the current
        % barycentre which is not in its class
        DistMat = DistMatAux(:,:,1).^0.5;
        
        [~, IidxNearestNeighbour] = sort(DistMat,1);
        
        IdxNextPointNotInClass = IidxNearestNeighbour(1, 1:nbarysInClass)';
            
        ClassNextPointNotInClass = ClassOfPointsAux(IdxNextPointNotInClass);
        PointInClass = BarycentresInClass(1:nbarysInClass,:);
        PointNotInClass = PointsNotInClass(IdxNextPointNotInClass,:);
        idxIclassVec = idxIclass*ones(nbarysInClass,1);

        [PointsLeftFinal, PointsRightFinal, Finished] = ...
            tripletsByBisection(PointInClass, PointNotInClass, ...
                                idxIclassVec, ClassNextPointNotInClass, ...
                                ProblemDescr, FaultApproxParams);

        for ipoint = 1: nbarysInClass
            % index of the point not in the current class, but nearest to
            % the given one
            if Finished(ipoint)
                NumPointsSurf{iclass, ClassValsInv(ClassNextPointNotInClass(ipoint))}(1) = ...
                    NumPointsSurf{iclass, ClassValsInv(ClassNextPointNotInClass(ipoint))}(1) + 1;

                NumPointsSurf{ClassValsInv(ClassNextPointNotInClass(ipoint)), iclass}(1) = ...
                    NumPointsSurf{ClassValsInv(ClassNextPointNotInClass(ipoint)), iclass}(1) + 1;

                iaux = ClassValsInv(ClassNextPointNotInClass(ipoint));
                PointSetsSurface{iclass, iaux}{1}(NumPointsSurf{iclass, iaux}(1),:) = ...
                    PointsLeftFinal(ipoint,:);
                PointSetsSurface{iaux, iclass}{1}(NumPointsSurf{iaux, iclass}(1),:) = ...
                    PointsRightFinal(ipoint,:);
            else
                warning('bisection failed in building block iniapprox.')
            end
        end
    end
    
    % Last part of iniapprox: remove duplicates and clusters. Two points
    % are considered duplicate if closer than FaultApproxParams.eps.
    % Moreover, remove almost duplicates: points that are closer to their
    % nearest neighbour than minDistFactor*maxDistForSurfacePoints are
    % collected into clusters wich are then reduced according to their
    % geometry (being in a cluster is transitive here).
    if ProblemDescr.verboseMode
        disp('-- remove duplicates in initial point sets')
    end
    for iclass = 1: nclasses -1

        for jclass = iclass + 1: nclasses
            [PointSetsSurface{iclass, jclass}{1}, PointSetsSurface{jclass, iclass}{1}] = ...
                removeDuplicates(PointSetsSurface{iclass, jclass}{1}, ...
                                 PointSetsSurface{jclass, iclass}{1}, ...
                                 FaultApproxParams.eps);
            
            NumPointsSurf{iclass, jclass} = size(PointSetsSurface{iclass, jclass}{1}, 1);
            NumPointsSurf{jclass, iclass} = size(PointSetsSurface{jclass, iclass}{1}, 1);
        end
    end


    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'initial';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end

    if ProblemDescr.verboseMode
        disp('-- remove clusters in initial point sets')
    end
    for iclass = 1: nclasses -1
        for jclass = iclass + 1: nclasses
    
            [PointSetsSurface{iclass,jclass}{1}, ...
             PointSetsSurface{jclass,iclass}{1}, ...
             NumPointsSurf{iclass,jclass}(1)] = ...
                removeClusters(PointSetsSurface{iclass,jclass}{1}, ...
                               PointSetsSurface{jclass,iclass}{1}, ...
                               NumPointsSurf{iclass,jclass}(1), ...
                               FaultApproxParams);

            NumPointsSurf{jclass,iclass}(1) = NumPointsSurf{iclass,jclass}(1);
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_cluster_rem';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end

    % Consistency check: if there are no points near any fault line at all
    % (for whatever reason), skip the computation.
    pointsOnFaultLines = false;
    for iclass = 1: nclasses -1
        for jclass = iclass + 1: nclasses
            pointsOnFaultLines = pointsOnFaultLines | ...
                                 any(PointSetsSurface{iclass, jclass}{1}(:,1));
        end
    end
    
    if (~pointsOnFaultLines)
        warning('Unable to find any points near a fault line. Skip computation.')
        
        return
    end
    
    if (ndim == 2)

        [PointSetsSurface, NumPointsSurf, ...
         LeftDomainStart, LeftDomainEnd, bsuccessful] = ...
            getTripletsNearFault2D(PointSetsSurface, NumPointsSurf, ...
                                   nclasses, ClassVals, ...
                                   FaultApproxParams, ProblemDescr);
        
    elseif (ndim == 3)

        [PointSetsSurface, NormalsSurface, ...
         NumPointsSurf, bsuccessful] = ...
            getTripletsNearFault3D(PointSetsSurface, NumPointsSurf, ...
                                   nclasses, ClassVals, ...
                                   FaultApproxParams, ProblemDescr);
    end
end