% Compute a rough approximation of fault lines based upon barycentres. We
% follow Allasia et al: Adaptive detection and approximation of unknown
% surface discontinuities from scattered data, Simulation in Modelling
% Practice and Theory, July 2009
% For a given sampling point, we consider the nNearestPoints nearest ones.
% If this subset contains points from more than one class, we compute the
% barycentres of the points in one class. As approximation to the fault
% line, we consider the arithmetic mean between all combinations of 
% barycentres within the nNearestPoints nearest points and store these
% new points in MeansOfBarycentres.
% In Allasia et al, the means of the barycentres are called barycentres
% as well. We deviate from this term.
%
% Input:
% - PointSet: Set of points (size npoints x ndim)
% - ClassOfPoints: classes of the points in PointSet. We assume that
%   the classes are ordered according to PointSet such that
%   ClassOfPoints(i) representes the class of the i-th point in PointSet
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - MeansOfBarycentres: (see description above)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function MeansOfBarycentres = getBarycentres(PointSet, ClassOfPoints, ...
                                             FaultApproxParameters)
    
    % number of points and dimension
    npoints = size(PointSet,1);
    ndim = size(PointSet,2);
    
    % number of nearest points to consider
    nNearestPoints = FaultApproxParameters.nNearestPoints; 
    
    % When removing duplicates, points closer than epsLoc are considered
    % identical.
    epsLoc = FaultApproxParameters.eps;
    
    % We do not know how many barycentres exist, but at most npoints ones.
    % We shorten the vector later on.
    MeansOfBarycentres = zeros(npoints, ndim);
    
    % total number of the new means of barycentres computed from PointSet
    nMeansOfBarycentres = 0;
    
    % lower nNearestPoints if necessary
    nNearestPoints = min(nNearestPoints, size(PointSet, 1));

    % Compute distance matrix for finding the nearest points.
    DistMatAux = zeros(npoints, npoints, ndim);
    for i = 1: npoints
        DistMatAux(:,i,:) = PointSet - PointSet(i,:);
    end
    
    DistMatAux = DistMatAux.^2;
    for idim = 2: ndim
        DistMatAux(:,:,1) = DistMatAux(:,:,1) + DistMatAux(:,:,idim);
    end
    DistMat = DistMatAux(:,:,1).^0.5;
    
    % matrix containing the indices of the nearest points for all points in
    % PointSet
    [~, IdxNearestNeighbours] = sort(DistMat,1);
    
    % loop over points
    for ipoint = 1: npoints
            
        % These are the nNearestPoints nearest points to the given point.
        IdxCurrentNearestNeighbours = IdxNearestNeighbours(1:nNearestPoints, ipoint);
        NearestNeighbours = PointSet(IdxCurrentNearestNeighbours,:);
            
        % A point is a fault point if in its neighborhood are points of
        % several classes.
        ClassOfNearestNeighbours = ClassOfPoints(IdxCurrentNearestNeighbours);
            
        % Find the values and the number of classes represented in
        % NearestNeighbours.
        ClassVals = unique(ClassOfNearestNeighbours);
        nclasses = size(ClassVals,1);
        
        % NearestNeighbours contains points from more than one class.
        if (nclasses > 1)

            % For each class, compute the barycentre of the points in
            % NearestNeighbours which belong to that class.
            BarycentresInClass = zeros(nclasses, ndim);
                            
            % Add barycentres for improved approximation of the fault.
            for iclass = 1: nclasses
                PointsInClass = NearestNeighbours(ClassOfNearestNeighbours == ClassVals(iclass), :);

                % arithmetic mean of the points within one class (aka
                % barycentres)
                BarycentresInClass (iclass, :)= sum(PointsInClass,1)/size(PointsInClass, 1);
            end
            
            % compute means of barycentres for any reasonable combination
            % of classes and store them in BarycentresInClass
            for iclass = 1: nclasses
                for jclass = iclass+1: nclasses
                    nMeansOfBarycentres = nMeansOfBarycentres + 1;
                    MeansOfBarycentres(nMeansOfBarycentres, :) = ...
                        0.5*(BarycentresInClass(iclass,:) + ...
                             BarycentresInClass(jclass,:));
                end
            end
                
        end
    end
    
    % Shorten vector to its appropriate length.
    MeansOfBarycentres = MeansOfBarycentres(1:nMeansOfBarycentres, :);

    % Last part: remove duplicates. Two points are considered duplicate if
    % closer than epsLoc.
    for i = 1: nMeansOfBarycentres-1
        for j = i+1: nMeansOfBarycentres
            if(abs(norm(MeansOfBarycentres(i,:) - MeansOfBarycentres(j,:))) < epsLoc)
                MeansOfBarycentres(j,:) = -42;
            end
        end
    end
    MeansOfBarycentres = MeansOfBarycentres(MeansOfBarycentres(:,1) ~= -42,:);
end