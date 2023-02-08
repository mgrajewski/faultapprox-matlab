% To remove clusters, it is important to remove the inner points, not
% the points on the boundary of the cluster. It is therefore
% straightforward in 2D to consider e.g. three consecutive points and
% to remove the one in the middle if these three points are very close.
% However, there is no sorting of the points along the fault line yet
% such that it is hard to tell which of the points is the one in the
% middle. Therefore, we perform local sorting. To do so, we assume that
% the fault line is somehow smooth such that we can find a new
% coordinate system in which the fault line locally resembles the
% x-axis. This coordinate system allows us then to establish sorting
% of the points by just comparing their (new) x-coordinates.
% We reduce the 3D-case to the 2D case. If there is a cluster, we have
% to find a representative point or, if the cluster is larger, some of
% them. Thanks to the eigendecomposition, we find the axis the cluster
% is oriented to and can sort the cluster points according to it. Then,
% we proceed as in 2D.
%
% Input:
% PointsIclass: n_points x n_dim-array of points in class iclass
% PointsJclass: n_points x n_dim-array of points in class jclass
% numPoints: number of points in PointsIclass (We assume that
%   PointIclass and PointJclass consist of the same number of points)
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% PointsIclass: n_points x n_dim-array of points in class iclass (after
%   cluster removal)
% PointsJclass: n_points x n_dim-array of points in class jclass (after
%   cluster removal)
% numPoints: number of points in PointsIclass (We assume that
%   PointIclass and PointJclass consist of the same number of points)
function [PointsIclass, PointsJclass, numPoints] = removeClusters(PointsIclass, ...
    PointsJclass, numPoints, FaultApproxParams)

    % dimension
    ndim = size(PointsIclass, 2);

    minDistFactor = FaultApproxParams.minDistFactor;
    maxDistForSurfacePoints = FaultApproxParams.maxDistForSurfacePoints;

    % compute distance matrix
    DistMatAux = zeros(numPoints, numPoints, ndim);
    for k = 1: numPoints
        DistMatAux(:,k,:) = PointsIclass - PointsIclass(k,:);
    end
    
    DistMatAux = DistMatAux.^2;
    for idim = 2: ndim
        DistMatAux(:,:,1) = DistMatAux(:,:,1) + DistMatAux(:,:,idim);
    end
    DistMat = DistMatAux(:,:,1).^0.5;

    % As the true number of points is decreased when removing clusters, we
    % store the original number of points.
    numPointsIni = numPoints;

    
    alreadyInCluster = zeros(1,numPointsIni);

    icluster = 0;
    ClusterPointsAll = {};
    for ipoint = 1: numPointsIni

        % point belongs to no cluster yet
        if (~alreadyInCluster(ipoint))
            ClusterPoints = DistMat(ipoint, :) < minDistFactor*maxDistForSurfacePoints;

            % The point itself has distance zero, but this is no cluster,
            % so there is no need to proceed.
            if (size(find(ClusterPoints),2) > 1)
                finished = false;
                while ~finished
                    ClusterPointsAux = zeros(1, numPointsIni);
                    iidx = find(ClusterPoints);
                    for jpoint = 2: size(iidx,2)
                        ClusterPointsNew = DistMat(iidx(jpoint), :) < minDistFactor*maxDistForSurfacePoints;
                        ClusterPointsAux = max(ClusterPointsAux, ClusterPointsNew);
                    end

                    finished = ~any(ClusterPoints - ClusterPointsAux);
                    ClusterPoints = max(ClusterPointsAux, ClusterPoints);
                end
                icluster = icluster + 1;
                alreadyInCluster = max(alreadyInCluster, ClusterPoints);
                ClusterPointsAll{icluster} = find(ClusterPoints);
            else
                alreadyInCluster(ipoint) = 1;
            end
        end    
    end
                        
    % remove clusters
    for icluster = 1: size(ClusterPointsAll,2)

        % If the cluster just consists of two points, take their mean
        % value.
        iidx = ClusterPointsAll{icluster};
        numPointsInCluster = size(iidx,2);
        if (numPointsInCluster == 2)

            MeanIclass = 1/numPointsInCluster*ones(1,numPointsInCluster)*PointsIclass(iidx,:);
            MeanJclass = 1/numPointsInCluster*ones(1,numPointsInCluster)*PointsJclass(iidx,:);
            PointsIclass(iidx(1),:) = MeanIclass;
            PointsJclass(iidx(1),:) = MeanJclass;
            PointsIclass(iidx(2),:) = -42;
            PointsJclass(iidx(2),:) = -42;

            numPoints = numPoints - 1;
        else

            % test, how large the cluster is
            DistMatLoc = DistMat(iidx, iidx);

            clusterSize = max(max(DistMatLoc));

            % the cluster is very small: take the mean
            if (clusterSize < minDistFactor*maxDistForSurfacePoints)
                MeanIclass = 1/numPointsInCluster*ones(1,numPointsInCluster)*PointsIclass(iidx,:);
                MeanJclass = 1/numPointsInCluster*ones(1,numPointsInCluster)*PointsJclass(iidx,:);
                PointsIclass(iidx(1),:) = MeanIclass;
                PointsJclass(iidx(1),:) = MeanJclass;
                PointsIclass(iidx(2:numPointsInCluster),:) = -42;
                PointsJclass(iidx(2:numPointsInCluster),:) = -42;

                numPoints = numPoints - numPointsInCluster + 1;

            % take selected points from the cluster: left most
            % and rightmost point and somewhere in the middle,
            % if necessary
            else
                numPointsMiddle = 0;

                % cluster points
                ClusterPoints = PointsIclass(iidx, :);

                % local midpoint
                xmid = 1/numPointsInCluster*ones(1, numPointsInCluster)*ClusterPoints;

                % compute local coordinate system based on the normal vector
                ClusterPoints = ClusterPoints - xmid;
                [~, ~, Q] = svd(ClusterPoints);
                
                % cluser points in the local coordinate system
                ClusterPoints = ClusterPoints*Q;
                
                % sort points locally according to the longest axis. As
                % singular values are provided in descending order, this is
                % alway the first coordinate.
                [~, IidxSorted] = sort(ClusterPoints(:,1));
                ClusterPointsSorted = ClusterPoints(IidxSorted, :);
                
                % cluster points in class iclass
                ClusterPointsIclass = PointsIclass(iidx(IidxSorted),:);

                % Replace the two cluster points with the lowest index
                % with the two extremal cluster points, as these will be
                % kept in any case.
                PointsIclass(iidx(1),:) = ClusterPointsIclass(1,:);
                PointsIclass(iidx(2),:) = ClusterPointsIclass(numPointsInCluster,:);

                ClusterPointsJclass = PointsJclass(iidx(IidxSorted),:);
                PointsJclass(iidx(1),:) = ClusterPointsJclass(1,:);
                PointsJclass(iidx(2),:) = ClusterPointsJclass(numPointsInCluster,:);

                % If the cluster is large, it may make sense to include
                % inner points as well. We start at one end and consider
                % the first point in the cluster which has sufficient
                % distance from the start point, if such a point exists.
                % If it is sufficiently remote from the end of the cluster,
                % we include this point. We proceed until there are no such
                % points left.
                stop = false;
                idxMid = [];
                istart = 1;
                while ~stop
                    % first point which is sufficiently far away from the
                    % first in the cluster (The points are sorted according
                    % to the x-coordinate after orthogonal transformation!)
                    idxPoint = find(ClusterPointsSorted(:,1) - ClusterPointsSorted(istart,1) > minDistFactor*maxDistForSurfacePoints, 1, 'first');

                    % there is no such point (may happen when some inner
                    % points have been found already)
                    if size(idxPoint, 1) == 0
                        stop = true;
                    else
                        % Test, if the first inner point sufficiently far
                        % away from the first point is sufficiently far
                        % away from the last point in the cluster. If so,
                        % it is a valid inner point.
                        if (ClusterPointsSorted(numPointsInCluster,1) - ClusterPointsSorted(idxPoint,1) < minDistFactor*maxDistForSurfacePoints)
                            stop = true;
                        else
                            idxMid = [idxMid, idxPoint];
                            ClusterPointsSorted(istart:idxPoint-1,1) = ClusterPointsSorted(idxPoint,1);
                            istart = idxPoint;
                            numPointsMiddle = numPointsMiddle + 1;
                        end
                    end
                end
                
                % mark all points except of the extremal points and
                % valid inner points for deletion
                PointsIclass(iidx(3+numPointsMiddle:numPointsInCluster),:) = -42;
                PointsIclass(iidx(3:3+numPointsMiddle-1),:) = ClusterPointsIclass(idxMid,:);

                PointsJclass(iidx(3+numPointsMiddle:numPointsInCluster),:) = -42;
                PointsJclass(iidx(3:3+numPointsMiddle-1),:) = ClusterPointsJclass(idxMid,:);

                numPoints = numPoints - numPointsInCluster + 2 + numPointsMiddle;
            end         
        end
    end
    PointsIclass = PointsIclass(PointsIclass(:,1) ~= -42,:);
    PointsIclass = PointsIclass(1:numPoints,:);
 
    PointsJclass = PointsJclass(PointsJclass(:,1) ~= -42,:);
    PointsJclass = PointsJclass(1:numPoints,:);
end