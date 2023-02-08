% The classification of points induces a subdivision of the domain
% into subdomains which themselves may consist of several components.
% faultApprox is to compute approximations to these sets.
% This subroutine returns the class and component index of the points
% in the array Points with respect to the approximations provided by
% faultApprox.
%
% This function uses the Matlab-library inpoly (see below) for the
% point-in-polygon tests. It is way faster than the Matlab-function
% inpolygon.
% Darren Engwirda (2020). INPOLY: A fast points-in-polygon test
% (https://www.github.com/dengwirda/inpoly), GitHub. 
%
% Input:
% - Points: (number of points x ndim)-array of point coordinates
% - Subdomains: Cell array containing polygonal approximations of the
%   subdomains and their components
%
% Output:
% - class: array (length: number of points) containing the classes of the
%   points in Points (-1, if the point could not be assigned to any class)
% - component: array (length: number of points) containing the component
%   indices of the classes the points in Points belong to (-1, if no
%   component could be found)

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [class, component] = classifyPoints(Points, Subdomains)
    ndim = size(Points, 2);
    if (ndim == 2)
        class = -1*ones(size(Points, 1),1);
        component = -1*ones(size(Points, 1),1);
        
        nclasses = size(Subdomains,1);
        for iclass = 1: nclasses
            for icomp = 1: size(Subdomains{iclass},1)
                in = inpoly2(Points, Subdomains{iclass}{icomp});
                class(in) = iclass;
                component(in) = icomp;
            end
        end
    end
end