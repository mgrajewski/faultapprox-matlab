% This function computes a triangle mesh suitable for use in Paraview from
% polygonal approximations to subdomains. It currently works in 2D only.
% The Mesh contains all points of the polygonal approximations of the
% subdomains twice. This allows for discontinuous representations in
% Paraview.
% The function relies for generating the 2D-mesh on Mesh2D, a library for
% grid generation written in Matlab by Darren Engwirda. The meshing
% algorithms are based upon
% D. Engwirda, (2014): "Locally-optimal Delaunay-refinement and
%   optimisation-based mesh generation", Ph.D. Thesis School of Mathematics
%   and Statistics, Univ. of Sydney.
%   http://hdl.handle.net/2123/13148
% It is available on https://github.com/dengwirda/mesh2d.
%
% Input:
% - Subdomains: Structure of polygonal approximations of subdomain
%   boundaries
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - ndim: dimension (2 or 3)
%
% Output:
% - GlobalPoints - All points of the mesh
% - GlobalMesh - Contains the global indices of the vertices of the
%   triangles
% - PointData - class of the points in GlobalPoints

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [GlobalPoints, GlobalMesh, PointData] = ...
    visualizeSubdivision(Subdomains, ProblemDescr, ndim)
                                                   
    global ncalls
                                                   
    % number of inner points in the subdomains
    % As mesh2d creates inner points on its own when necessary, we can set
    % it to 0 as long we use mesh2D. However, if we use the Matlab Delaunay
    % mesher, it should be set to some positive value.  
    nInnerPoints = 0;
    
    if (ndim == 0)
        error('The Subdomains-structure appears to be empty.')
    elseif (ndim == 2)
        % visualisation
        SubdomainMeshConnect = cell(size(Subdomains,1),1);
        SubdomainMeshVerts = cell(size(Subdomains,1),1);

        nverts = 0;
        ncells = 0;
        ncallsOld = ncalls;

        for iclass = 1: size(Subdomains,1)

            for icomp = 1 : size(Subdomains{iclass},1)

                % create some auxiliary points inside the subdomains to ease
                % triangulation
                
                % remove all NaNs
                bdryPoints = Subdomains{iclass}{icomp}(~isnan(Subdomains{iclass}{icomp}(:,1)),:);

                % number of actual boundary points
                npoints = size(bdryPoints,1);
                
                % mean value
                mean = bdryPoints'*ones(npoints,1);
                mean = mean'./npoints;
                
                % construct bounding box; orientation by SVD
                [~,~,V] = svd(bdryPoints-kron(mean, ones(npoints,1)));
                PointsAux = bdryPoints*V;
                sizes = max(PointsAux) - min(PointsAux);
                mean = min(PointsAux) + 0.5*sizes;
                mean = mean*V';
                
                % this is on [-0.5, 0.5]^ndim
                PointSetVis = CreateHaltonSet(nInnerPoints, ndim) - 0.5;

                % transform the points according to the bounding box
                % we want points inside, therefore the factor 0.85
                PointSetVis = PointSetVis*0.85*diag(sizes)*V' + kron(mean, ones(nInnerPoints,1));

                % test, if these points are really inside the subdomain
                ClassOfPointsVis = computeClassification(PointSetVis, ProblemDescr);

                % MATLAB is one-based, however, in paraview, we need it 0-based
                aux = size(Subdomains{iclass}{icomp},1)-1;

                % if the components are not simply connected, there are several
                % boundary components. Following a MATLAB convention, they are
                % separated by NaNs in the node vector. Unfortunately, Mesh2D
                % does not follow that convention such we have to remove all
                % NaNs and build the edge vectors accordingly.
                IdxNaNs = find(isnan(Subdomains{iclass}{icomp}(:,1)));
                IdxNaNs = [0 IdxNaNs aux+1];
                edges = [];
                for i = 1: size(IdxNaNs,2)-1
                    edges = [edges; (IdxNaNs(i)+1:IdxNaNs(i+1)-2)'-i+1, (IdxNaNs(i)+1:IdxNaNs(i+1)-2)'-i+2];
                    edges = [edges; IdxNaNs(i+1) - i, IdxNaNs(i)+2-i];
                end

                % get all points in the current class
                AuxPoints = PointSetVis(ClassOfPointsVis == iclass, :);

                % take the ones which are inside the polygon approximating
                % the boundary of the subdomain
                in = inpolygon(AuxPoints(:,1), AuxPoints(:,2), ...
                               Subdomains{iclass}{icomp}(:,1), ...
                               Subdomains{iclass}{icomp}(:,2));
                AuxPoints = AuxPoints(in,:);

                % remove the NaNs and add the inner points
                SubdomainMeshVerts{iclass}{icomp} = Subdomains{iclass}{icomp}(~isnan(Subdomains{iclass}{icomp}(:,1)),:);
                SubdomainMeshVerts{iclass}{icomp} = [SubdomainMeshVerts{iclass}{icomp}(1:end-1,:); AuxPoints];

                % Workaround should be obsolete for mesh2D 3.2.2
                % call external library for meshing (workaround for a bug in
                % mesh2D which cannot deal with a simple triangle as input
    %            if (size(SubdomainMeshVerts{iclass}{icomp},1) > 3)
                    [SubdomainMeshVerts{iclass}{icomp}, ~, ...
                     SubdomainMeshConnect{iclass}{icomp}, ~] = ...
                     refine2(SubdomainMeshVerts{iclass}{icomp}, edges);
    %            else
    %                SubdomainMeshConnect{iclass}{icomp} = [1 2 3];
    %            end

                nverts = nverts + size(SubdomainMeshVerts{iclass}{icomp},1);
                ncells = ncells + size(SubdomainMeshConnect{iclass}{icomp},1);
            end
        end

        ncalls = ncallsOld;

        GlobalMesh = zeros(ncells,4);
        GlobalPoints = zeros(nverts,2);
        PointData = zeros(nverts,1);
        istartConnect= 1;
        istartVerts = 1;
        for iclass = 1:size(Subdomains,1)
            for icomp = 1 : size(Subdomains{iclass},1)
                iendConnect = istartConnect + size(SubdomainMeshConnect{iclass}{icomp},1)-1;
                iendVerts = istartVerts + size(SubdomainMeshVerts{iclass}{icomp},1)-1;
                PointData(istartVerts:iendVerts) = iclass;

                GlobalMesh(istartConnect:iendConnect,1:3) = SubdomainMeshConnect{iclass}{icomp} + istartVerts-1;
                GlobalPoints(istartVerts:iendVerts,:) = SubdomainMeshVerts{iclass}{icomp};
                istartConnect = iendConnect + 1;
                istartVerts = iendVerts + 1;
            end
        end

        % ParaView is 0-based, whereas the indices are (up to here) 1-based!
        GlobalMesh(:, 1:3) = GlobalMesh(:, 1:3)-1;
        GlobalMesh(:,4) = 3;
    elseif (ndim == 3)
        ncells = size(Subdomains{1}{2}{1},1);
        nverts = size(Subdomains{1}{2}{2},1);
        GlobalMesh = zeros(ncells, 4);
        GlobalMesh(:, 1:3) = Subdomains{1}{2}{1} - 1;
        GlobalMesh(:,4) = 3;
        GlobalPoints = Subdomains{1}{2}{2};
        PointData = ones(nverts, 1);
    else
        error(['The dimension of the data is ',int2str(ndim), ...
               '. This function supports only dimensions 2 and 3.'])
    end
end