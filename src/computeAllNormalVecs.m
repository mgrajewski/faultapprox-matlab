% This function computes the approximate normals for th set of triplets
% given by PointsIclass and PointsJclass, respectively. The normals are
% based on local linearization in compueNormalVec and point to jclass.
%
% Input:
% - PointsIclass, PointsJclass: set of triplets
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - NormalsSurface: approximate normal vector for each triplet

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function NormalsSurface = computeAllNormalVecs(PointsIclass, PointsJclass, ...
                                               FaultApproxParams)

    [npoints, ndim] = size(PointsIclass);
    nNearestPoints = min(npoints, FaultApproxParams.nNearestPoints);

    % vector to contain all normal vectors
    NormalsSurface = zeros(npoints,3);
    
    if (npoints >= ndim)

        % compute the distance matrix (actually the squares of
        % distances)
        DistMat = zeros(npoints, npoints, ndim);
        for k = 1: npoints
            DistMat(:,k,:) = PointsIclass-PointsIclass(k,:);
        end

        DistMat = DistMat.^2;
        for idim = 2: ndim
            DistMat(:,:,1) = DistMat(:,:,1) + DistMat(:,:,idim);
        end
        [~, IidxNearestNeighbour] = sort(DistMat(:,:,1),1);

        % We consider the nNearestPoints nearest neighbours and utilize
        % them to approximate a normal vector using computeNormalVec.
        for ipoint = 1: npoints
            
            NearestPoints = PointsIclass(IidxNearestNeighbour(1:nNearestPoints, ipoint),:);
            NormalsSurface(ipoint,:) = computeNormalVec(NearestPoints);

            % The difference PointsJclass(ipoint,:) -
            % PointsIclass(ipoint,:) is a crude approximation to the true
            % normal, but has the correct orientation in any case.
            % Therefore, the angle between that approximation and the true
            % normal is 90° at most, indicated by a non-negative scalar
            % product. If so, the normal vector has the wrong orientation
            % and has to be flipped.
            if(NormalsSurface(ipoint,:)*(PointsJclass(ipoint,:) - PointsIclass(ipoint,:))' < 0)
                NormalsSurface(ipoint, :) = -NormalsSurface(ipoint,:);
            end
        end
    else
        warning(['Computing normals requires the point set to consist ' ...
                 'of at least ' int2str(ndim) ' points in ' ...
                 int2str(ndim) ' dimensions.'])
    end
end