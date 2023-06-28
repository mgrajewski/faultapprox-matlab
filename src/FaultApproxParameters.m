% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
classdef FaultApproxParameters
        
    properties

        % points closer as eps are regarded as identical
        eps = 1e-10;
        
        % number of the nearest points to consider when computing
        % barycentres
        nNearestPoints = 10;
        
        % GetPointsNearSurface
        
        % desired maximum distance of a point on the fault line to the next one
        maxDistForSurfacePoints = 0.05;

        % for safeguarding in computing valid point pairs as starting
        % values for bisection
        alpha = 0.25;

        % points closer than minDistFactor*maxDistForSurfacePoints are removed
        % if appropriate
        minDistFactor = 0.2;

        % cosine of maximal admissible angle between consecutive line segments for sorting
        cosAlphaMax = -0.9;
        
        % number of nearest points to consider for sorting
        nNearestSorting = 5;
                
        % maximal number of adaptive refinement/coarsening steps for
        % adaptively placing points on the fault line/surface
        maxiterAdapt = 4;
        
        % number of points for computing a local coordinate system and
        % normal vectors in case of a fault line
        NumPointsLocal = 3;
                
        % line segments with error lower than errMin are coarsened
        errMin = 0.001;
        
        % line segments with error larger than errMax are refined
        errMax = 0.15;
        
        % parameters for singleTripletByBisection
        
        %stopping criterion for bisection algorithm
        abstolBisection = 1e-3;

        % maximum number of bisection iterations
        maxiterBisection = 60;
        
        % There may be a large distance between neighbouring points
        % near a fault line. We try to fill such gaps by adding points.
        % However, to the additional information provided, the (assumed)
        % shape of the fault line may changes considerably. Due to that,
        % new gaps may emerge. Therefore, we put the filling process in a
        % loop and perform at most maxTrialsForFillingGapsInLines of such
        % steps.
        %
        % In 3D, it takes sometimes several passes to really fill all holes
        % in the representation of a surface. This is the number of passes.
        maxTrialsForFillingGaps = 3;
        
    end
end