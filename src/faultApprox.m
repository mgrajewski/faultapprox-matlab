% We consider piecewise smooth a function f: \Omega -> {1,..., n} such that
% there are lines or more generally manifolds where the function is
% discontinuous, and \Omega is an axis-parallel rectangle or a cuboid.
% We assume that for \Omega_i = f^{-1}(i),
%   \Omega = \overline{\Omega_i} \cup \hdots \overline{\Omega_n}.
% The purpose of this function is to represent such lines/surfaces by
% sufficiently many points in their near vicinity. We will call a curve
% or surface \overline{\Omega_i} \cap \overline{\Omega_j} a fault in what
% follows.
% We proceed as follows:
% 1) Starting from an initial sample set, we detect how many of such
%    subdomains exist.
% 2) step initialise: We detect all sample points in the vicinity of a
%    fault surface and utilise them to find additional points nearer to
%    the faults based on barycentres. We repeat that process once. We
%    classify these points accordingly and use for them for finding
%    additional ones extremely close to the decision line/surface with a 
%    bisection approach (Building block iniapprox).
% 3) step fill: We detect gaps or holes in the representing point set,
%    which we fill in this step.
% 4) step expand: We ensure that the set of points represents the fault in
%    its entire extent.
% 5) step adapt: We adaptively add or remove (2D only) points in the
%    representing sets depending on geometric properties like curvature. We
%    end up with a point cloud representing the faults sufficiently
%    accurately.
% 6) step reconstruct: In 2D, we describe each of the subdomains \Omega_i
%    by a closed polygonal approximation of its boundary based on the
%    corresponding faults. The points the polygon consists of are ordered
%    counter-clockwise such that the subdomain appears on the right side
%    when following the points.
% 7) We provide a visualisation as vtu-file if desired.
%
% Input:
% - PointSet: initial set of points where to sample f; Dimension is
%   number of points  x ndim, where ndim is the dimension of \Omega, 2 or
%   3.
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - Subdomains: Cell array of arrays. Its entries contain the polygons
%   of the subdomains as y x ndim-arrays. As the polygons are closed,
%   the first and last vertex coincide.
% - reconstructionSucceeded: boolean. True, if (in 2D), a closed polygonal
%   approximation of any \Omega_i could be obtained.
% - ncallsResults: number function evaluations necessary in the algorithm;
%   ncallsResult(1): number of callings of the function
%   ncallsResult(2): number of evaluations of f
% - tel: elapsed computational time

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [Subdomains, reconstructionSucceeded, ncallsResult, tel] = ...
    faultApprox(PointSet, ProblemDescr, FaultApproxParams)

    tstart = tic;

    % Using a global variable for counting ncalls is not up to high
    % standards, but makes the parameter lists less cluttered.
    global ncalls;

    if ProblemDescr.extendedStats
        global ExtendedStats;
        if ~isa(ExtendedStats, 'Statistics')
            ExtendedStats = Statistics;
        end
    end

    ncalls = [0 0];
    ndim = size(PointSet, 2);
    
    % Remove points which are obviously not inside the domain to consider.
    for idim = 1: ndim
        PointSet = PointSet(PointSet(:, idim)<= ProblemDescr.Xmax(idim),:);
        PointSet = PointSet(PointSet(:, idim)>= ProblemDescr.Xmin(idim),:);
    end

    % Classification of the initial points.
    ClassOfPoints = computeClassification(PointSet, ProblemDescr);

    % ClassVals contains the class values found so far. These are not
    % necessarily the class indices: Imagine that f(\Omega) = {1,2,5} and
    % all corresponding subdomains have been covered by PointSet. Then,
    % ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
    ClassVals = unique(ClassOfPoints);

    % If there is more than one class present in the initial point set, we
    % start reconstructing the subdomains.
    if size(ClassVals,1) > 1

        if ProblemDescr.verboseMode
            disp('- compute set of barycentres M')
        end

        % First rough approximation of the faults by means of barycentres.
        MeansOfBarycentres = getBarycentres(PointSet, ClassOfPoints, ...
                                            FaultApproxParams);
        
        % classes of the means of barycentres
        ClassMeansOfBarycentres = computeClassification(MeansOfBarycentres, ...
                                                        ProblemDescr);
        
        if ProblemDescr.extendedStats
            ExtendedStats.pos_in_code{end+1} = 'after_M';
            ExtendedStats.MeansBarys{end+1} = MeansOfBarycentres;
            ExtendedStats.ClassMeansBarys{end+1} = ClassMeansOfBarycentres;
            ExtendedStats.ncalls{end+1} = ncalls(2);
        end
        % Second slightly less rough approximation of the faults by means
        % of barycentres of means of barycentres.
        if ProblemDescr.verboseMode
            disp('- compute second set of barycentres M2')
        end
        MeansOfBarycentres2 = getBarycentres(MeansOfBarycentres, ...
                                             ClassMeansOfBarycentres, ...
                                             FaultApproxParams);
        ClassMeansOfBarycentres2 = computeClassification(MeansOfBarycentres2, ...
                                                         ProblemDescr);

        if ProblemDescr.extendedStats
            ExtendedStats.pos_in_code{end+1} = 'after_M2';
            ExtendedStats.MeansBarys2{end+1} = MeansOfBarycentres2;
            ExtendedStats.ClassMeansBarys2{end+1} = ClassMeansOfBarycentres2;
            ExtendedStats.ncalls{end+1} = ncalls(2);
        end

        % It may happen that not all classes are detected in the initial
        % point set, but only later in the extended point set classified so
        % far. Therefore, we recompute ClassVals here.
        ClassVals = unique([ClassOfPoints; ClassMeansOfBarycentres; ...
                            ClassMeansOfBarycentres2]);
        
        % getTripletsNearFaults contains the steps iniapprox, fill, expand
        % and adapt.
        if ProblemDescr.verboseMode
            disp('- compute points on subdomain boundaries')
        end
        [PointSetsSurface, LeftDomainStart, LeftDomainEnd, ...
         bsuccessful, NormalsSurface] = ...
            getTripletsNearFaults([MeansOfBarycentres; MeansOfBarycentres2], ...
                                  [PointSet; MeansOfBarycentres; MeansOfBarycentres2], ...
                                  [ClassMeansOfBarycentres; ClassMeansOfBarycentres2], ...
                                  [ClassOfPoints; ClassMeansOfBarycentres; ...
                                   ClassMeansOfBarycentres2], ...
                                  ProblemDescr, FaultApproxParams);

        if (~bsuccessful)
            reconstructionSucceeded = false;
            Subdomains = 0;
            ncallsResult = [0 0];
            tel = 0;
            return;
        end        
        
        % Step reconstruct: Describe the subdomains approximately as closed
        % polygons according to the classes using the points near the
        % faults (only for 2D).
        if ProblemDescr.verboseMode
            disp('- reconstruct subdomains')
        end
        [Subdomains, reconstructionSucceeded] = ...
            reconstructSubdomains(PointSetsSurface, LeftDomainStart, ...
                                  LeftDomainEnd, ProblemDescr, ...
                                  ClassVals, FaultApproxParams);
        tel = toc(tstart);

        ndim = size(PointSet,2);
        % In case of three dimensions: abuse the Subdomains-Structure for
        % returning the outer normal vector for each point.
        if (ndim == 3)
            for iclass = 1: size(ClassVals,1)
                for jclass = iclass+1: size(ClassVals,1)

                    Subdomains{iclass}{jclass}{3} = NormalsSurface{iclass,jclass}{1};
                end
            end
        end
        
        % Write the result in a vtu-file for visualisation, if desired.
        if (~strcmp(ProblemDescr.OutputFileVTU ,'') && reconstructionSucceeded)

            % prepare visualisation
            [GlobalPoints, GlobalMesh, PointData] = ...
                visualizeSubdomains(Subdomains, ProblemDescr, ndim);
           
            if ProblemDescr.verboseMode
                disp('- write visualisation file')
            end
            % export as vtu-file for further analysis in e.g. paraview
            Export2VTU(GlobalPoints, GlobalMesh, PointData, {'class'}, ...
                                       ProblemDescr.OutputFileVTU)
        end

    % All points in the initial point set belong to the same class: display
    % a warning and return dummy values.
    else
        warning('All points in the initial point set belong to the same class.')                       
        warning('Consider employing an enriched set of sampling points.')

        Subdomains = cell(ClassVals,1);
        aux = ProblemDescr.Xmin(1) - 1;
        epsLoc = 1e-9;
        for iclass = 1: ClassVals
            Subdomains{iclass} = cell(1);

            % dummy domain outside the real domain
            Subdomains{iclass}{1} = [aux aux; aux+epsLoc aux; ...
                                     aux+epsLoc aux+epsLoc; aux aux]; 
        end
        Subdomains{ClassVals}{1} = [ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmin(2)];
        reconstructionSucceeded = true;
        tel = toc(tstart);
    end
    ncallsResult = ncalls;
end

