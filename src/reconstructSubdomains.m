% The purpose of this function is to combine components of fault lines to
% polygons which approximate the subdomains associated with classification.
% 
% Input:
% - PointSetsSurface: Data structure containing the points near the fault
%   lines. PointSetsSurface{iclass, jclass} points to a structure of
%   arrays (as many as the components the fault line consists of)
%   containing the coordinates of the points near the fault
%   line between classes iclass and jclass which are themselves belong to
%   iclass.
% - LeftDomainStart: Data structure containing the index of the domain
%   boundary from which a component of a fault line starts (0, if it
%   starts in the interior).
%   indices of the domain edges:
%             3
%      _______________ 
%     |               |
%     |               |
%     |               |
%   4 |               | 2
%     |               |
%     |               |
%     |_______________|
%             1
%    
% - LeftDomainEnd: same as LeftDomainStart, but referring to the domain
%   boundary where a component of a fault line ends
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - ClassVals: Array containing the class values. These are not
%   necessarily the class indices. Imagine that f(\Omega) = {1,2,5}. Then,
%   ClassVals = [1,2,5], whereas the class indices range from 1 to 3.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% Subdomains: Structure polygonal approximation to the components of the
% subdomains associated to the several classes    

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [Subdomains, succeeded] = reconstructSubdomains(PointSetsSurface, LeftDomainStart, ...
                                                         LeftDomainEnd, ProblemDescr, ...
                                                         ClassVals, FaultApproxParams)
    
    
    % number of occurring classes
    nclasses = size(PointSetsSurface, 1);

    % find out dimension
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses
            PointSet = PointSetsSurface{iclass, jclass};
            if (~isempty(PointSet))
                dim = size(PointSet{1},2);
                break;
            end
        end
    end
    
    maxClassIdx = max(ClassVals);

    % create Subdomains-structure
    Subdomains = cell(maxClassIdx,1);

    if (dim == 2)

        succeeded = false;

        % maximal number of trials for completing a subdomain
        maxTrials = 100;


        PointSetsSurfaceAux = cell(nclasses);

        LeftDomainStartAux = -ones(size(LeftDomainStart));
        LeftDomainEndAux = -ones(size(LeftDomainEnd));

        % PointSetsSurface{iclass, jclass} and
        % PointSetsSurface{jclass, iclass} both contain points near the
        % fault line between the classes iclass and jclass. The points of 
        % the former one belong to iclass, the points of the latter one to
        % jclass. While this is worthwile for computing approximations of 
        % the components of a fault line, it is not for providing a
        % consistent subdivision of the domain by polygons which we intend
        % here. Therefore, we take the mean of these pairs of points.
        for iclass = 1:nclasses
            for jclass = 1:iclass-1
                ncomps = size(PointSetsSurface{iclass, jclass},1);
                for icomp = 1 : ncomps
                    PointSetsSurfaceAux{iclass, jclass}{icomp} = ...
                        0.5*(PointSetsSurface{iclass, jclass}{icomp} + flip(PointSetsSurface{jclass, iclass}{icomp}));
                    LeftDomainStartAux(iclass, jclass, icomp) = LeftDomainEnd(jclass, iclass, icomp);
                    LeftDomainEndAux(iclass, jclass, icomp) = LeftDomainStart(jclass, iclass, icomp);
                end
            end
        end
        
        for iclass = 1:nclasses
            for jclass = iclass+1:nclasses
                ncomps = size(PointSetsSurface{iclass, jclass},1);
                for icomp = 1 : ncomps
                    PointSetsSurfaceAux{iclass, jclass}{icomp} = flip(PointSetsSurfaceAux{jclass, iclass}{icomp});
                    LeftDomainStartAux(iclass, jclass,icomp) = LeftDomainStart(iclass, jclass, icomp);
                    LeftDomainEndAux(iclass, jclass, icomp) = LeftDomainEnd(iclass, jclass, icomp);
                end
            end
        end

        
        % Consistency check: if a component of a subdomain boundary starts
        % inside the domain but no other component related to this
        % subdomain end inside the domain, at least one boundary component
        % must be missing.
        for iclass = 1: nclasses

            % There are components belonging to iclass ending inside the
            % domain.
            iend = LeftDomainEndAux(iclass,:, :);
            istart = LeftDomainStartAux(iclass,:, :);

            % Something ends inside but nothing starts inside or something
            % starts inside but nothing ends inside.
            if (any(iend(:) == 0) && ~any(istart(:) == 0)) || (~any(iend(:) == 0) && any(istart(:) == 0))
                warnMsg = ['Some components of the boundary of class ' int2str(ClassVals(iclass)) 'are missing.'];
                warning(warnMsg);

                warnMsg = 'Reconstruction of subdomains failed. Consider taking a finer initial point set.';
                warning(warnMsg);
                return;
            end
        end


        % compute the maximal number of components per fault line
        maxNumComps = 0;
        for iclass = 1: nclasses
            for jclass = iclass + 1: nclasses
                maxNumComps = max(maxNumComps, size(PointSetsSurfaceAux{iclass, jclass},2));
            end
            Subdomains{ClassVals(iclass)} = cell(maxNumComps,1);
        end

        ParValsStart = -ones(nclasses, nclasses, maxNumComps);
        ParValsEnd = -ones(nclasses, nclasses, maxNumComps);

        % Allocate and fill NumPointsSurf with the number of points per
        % component. As the points in (iclass, jclass) are the same as in
        % (jclass, iclass) by construction, we exploit symmetry.
        NumPointsSurf = zeros([ size(PointSetsSurfaceAux) maxNumComps]);
        for iclass = 1: nclasses
            for jclass = iclass+1: nclasses
                for icomp = 1: size(PointSetsSurfaceAux{iclass, jclass},2)
                    NumPointsSurf(iclass, jclass, icomp) = size(PointSetsSurfaceAux{iclass, jclass}{icomp},1);
                    NumPointsSurf(jclass, iclass, icomp) = size(PointSetsSurfaceAux{jclass, iclass}{icomp},1);
                end
            end
        end
            %    direction of parametrisation of the global domain
        %                <-----------------
        %   |-----------------------------------------|
        %   |            <-----               <-----  |
        %   |       par subdomain I           ---->   |
        %   |                           _   /---------| X
        %   |                           /| /  <----   |
        %   |                          /  / par II    |
        %   |                         /  /         /\ |    /\
        %   |              I            |          |  |    |
        %   |                           |    II,1  |  |    |
        %   |                            \            |    |
        %   |                              \   ---->  |    |
        %   |                                \________| Y
        %   |                                         |
        %   |                                   /-----| Z
        %   |                                  /      |
        %   |                                 / II,2  |
        %   |       ----->                   /________|
        %   |_________________________________________|
        %
        %
        % The surface point sets for (iclass, jclass) are ordered such that
        % subdomain iclass is right when running along the fault line. This,
        % however is not sufficient for assembling the subdomains, if some
        % of them consist of several components: We start assembling
        % subdomain I by starting in X and end with this component of the
        % fault line in Y. Now, we need to continue with the second
        % component of the fault line starting in Z. To do so, we search
        % for fault line components starting at the right domain boundary
        % and find the one starting in Z (if there are more than one
        % starting at the right domain boundary, we take the closest one).
        % Realising that there are no more components, we close the
        % subdomain and are finished.
        % Now, we assemble both components of subdomain II. We may start in
        % Y and end in X. In contrast to above, we must not consider any
        % additional boundary components but instead stop and close that
        % first component of subdomain II. Ths same holds for component
        % II,2. 
        %
        % The main difference between these two cases is that for subdomain
        % II, the parameter value of the starting point is lower than the
        % parameter value of the end point, both being boundary points.
        % This is why we need parameter values for start and end points of
        % boundary components if they are boundary points.
        for iclass = 1: nclasses
            for jclass = 1: nclasses
                for icomp = 1: size(PointSetsSurfaceAux{iclass, jclass},2)
                    switch LeftDomainStartAux(iclass, jclass, icomp)
                        case(1)
                            aux = 1/(ProblemDescr.Xmax(1) - ProblemDescr.Xmin(1));
                            ParValsStart(iclass, jclass, icomp) = (PointSetsSurfaceAux{iclass, jclass}{icomp}(1,1) - ProblemDescr.Xmin(1))*aux;                      
                        case(2)
                            aux = 1/(ProblemDescr.Xmax(2) - ProblemDescr.Xmin(2));
                            ParValsStart(iclass, jclass, icomp) = 1 + (PointSetsSurfaceAux{iclass, jclass}{icomp}(1,2) - ProblemDescr.Xmin(2))*aux;                        
                        case(3)
                            aux = 1/(ProblemDescr.Xmax(1) - ProblemDescr.Xmin(1));
                            ParValsStart(iclass, jclass, icomp) = 3 - (PointSetsSurfaceAux{iclass, jclass}{icomp}(1,1)- ProblemDescr.Xmin(1))*aux;                        
                        case(4)
                            aux = 1/(ProblemDescr.Xmax(2) - ProblemDescr.Xmin(2));
                            ParValsStart(iclass, jclass, icomp) = 4 - (PointSetsSurfaceAux{iclass, jclass}{icomp}(1,2) - ProblemDescr.Xmin(2))*aux;
                    end

                    switch LeftDomainEndAux(iclass, jclass, icomp)
                        case(1)
                            aux = 1/(ProblemDescr.Xmax(1) - ProblemDescr.Xmin(1));
                            ParValsEnd(iclass, jclass, icomp) = (PointSetsSurfaceAux{iclass, jclass}{icomp}(NumPointsSurf(iclass, jclass, icomp))- ProblemDescr.Xmin(1))*aux;             
                        case(2)
                            aux = 1/(ProblemDescr.Xmax(2) - ProblemDescr.Xmin(2));
                            ParValsEnd(iclass, jclass, icomp) = 1 + (PointSetsSurfaceAux{iclass, jclass}{icomp}(NumPointsSurf(iclass, jclass, icomp),2) - ProblemDescr.Xmin(2))*aux;                        
                        case(3)
                            aux = 1/(ProblemDescr.Xmax(1) - ProblemDescr.Xmin(1));
                            ParValsEnd(iclass, jclass, icomp) = 3 - (PointSetsSurfaceAux{iclass, jclass}{icomp}(NumPointsSurf(iclass, jclass, icomp),1)) - ProblemDescr.Xmin(1)*aux;                        
                        case(4)
                            aux = 1/(ProblemDescr.Xmax(2) - ProblemDescr.Xmin(2));
                            ParValsEnd(iclass, jclass, icomp) = 4 - (PointSetsSurfaceAux{iclass, jclass}{icomp}(NumPointsSurf(iclass, jclass, icomp),2) - ProblemDescr.Xmin(2))*aux;
                    end             
                end
            end
        end

        CornersNext = zeros(4,2);
        CornersNext(1,:) = [ProblemDescr.Xmax(1) ProblemDescr.Xmin(2)];
        CornersNext(2,:) = [ProblemDescr.Xmax(1) ProblemDescr.Xmax(2)];
        CornersNext(3,:) = [ProblemDescr.Xmin(1) ProblemDescr.Xmax(2)];
        CornersNext(4,:) = [ProblemDescr.Xmin(1) ProblemDescr.Xmin(2)];


        % Combine components of fault line to polygons which approximate
        % the components of the subdomain associated to the classes.
        for iclass = 1: nclasses

            % We do not know how many components the subdomain associated
            % with iclass consists of, but at most as many as components as
            % the fault lines consist of.
            for icomp = 1 : maxNumComps

                % true, if all components of a subdomain are found
                bfinishedSub = false;

                % true, if all components of a component are found
                bfinishedPart = false;

                binitial = true;

                parStart = -1;
                ParEnd = 0;

                % If there are any components of a fault line associated
                % with iclass (we abuse NumPointsSurf to indicate that).
                if (any(any(NumPointsSurf(iclass, :,:) > 0)))

                    % Find the indices of the classes on the other side
                    % jclass and the component index jcomp of all these
                    % fault line components (may be several ones).
                    AuxArr = squeeze(NumPointsSurf(iclass, :, :) > 0);
                    if (size(AuxArr,1) == 1)
                        AuxArr = AuxArr';
                    end
                    [jclass, jcomp] = ind2sub(size(AuxArr), find(AuxArr > 0, 1));
                else
                    warning(['There are no known points on the boundary of class ' num2str(ClassVals(iclass))]);
                    succeeded = false;
                    return
                end

                % to avoid that this part is taken once more later on
                LeftDomainStartAuxInit = LeftDomainStartAux(iclass, jclass, jcomp);
                itry = 1;
                while ~bfinishedPart
                    itry = itry+1;

                    if (itry > maxTrials)
                        warning("Reconstruction of subdomain " + int2str(ClassVals(iclass)) + " failed.")
                        succeeded = false;
                        return
                    end

                    if jclass > 0
                        % it is not entirely clear why this error case
                        % should occur, but it does.
                        if (size(Subdomains{ClassVals(iclass)},1) >= icomp)
                            Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; PointSetsSurfaceAux{iclass, jclass}{jcomp}];
                        else
                            warning("Reconstruction of subdomain " + int2str(ClassVals(iclass)) + " failed.")
                            succeeded = false;
                            return
                        end                            
                        if binitial
                            edgeIni = LeftDomainStartAux(iclass, jclass, jcomp);
                            binitial = false;
                            parStart = ParValsStart(iclass, jclass, jcomp);
                        end
                        LeftDomainStartAux(iclass, jclass, jcomp) = -42;
                        edgeEnd = LeftDomainEndAux(iclass, jclass, jcomp);
                        NumPointsSurf(iclass, jclass, jcomp) = 0;
                    end

                    % If there are no appropriate boundary parts at all,
                    % stop reconstructing this subdomain and just close it.
                    if ~any(NumPointsSurf(iclass, :) )
                        bfinishedSub = true;
                        bfinishedPart = true;

                    % There are still unassigned boundary components
                    % belonging to the current subdomain.
                    else

                        jclass = 0;

                        % The last assigned boundary part ends on the domain
                        % boundary part with index edgeEnd.
                        if (edgeEnd > 0)

                            % Test if any other still unassigned boundary
                            % part starts on the current domain boundary.
                            % squeeze transforms the 1 x m x n-3D-array to
                            % an ordinary matrix.
                            iidx = squeeze(LeftDomainStartAux(iclass, :, :) == edgeEnd);
                            if (size(iidx,1) == 1)
                                iidx = iidx';
                            end

                            % Get x- or y-coordinate from the end point of
                            % the current boundary part (the one which ends
                            % on the current edge). This point is the last
                            % in the array of boundary points.
                            coordLast = Subdomains{ClassVals(iclass)}{icomp}(end,mod(edgeEnd+1,2)+1);
                            parLast = compParValues(edgeEnd, coordLast);

                            % There might be several boundary parts of that
                            % kind. Take the "nearest one" aka the one with the
                            % lowest starting parameter but still larger than the
                            % present end parameter.
                            % If there is no such boundary part, take the
                            % appropriate part of the domain boundary. Test
                            % before if the component of the subdomain is
                            % closed yet.
                            if (any(iidx(:)))
                                 % get the indices of the corresponding boundary
                                % parts
                                [i,j] = ind2sub(size(iidx), find(iidx ==1));
                                iidx2 = zeros(size(i,1), 2);
                                iidx2(:,1) = i;
                                iidx2(:,2) = j;

                                auxVals = zeros(size(iidx2, 1),1);
                                 % find the nearest part
                                for i = 1: size(iidx2, 1)
                                    auxVals(i) = PointSetsSurfaceAux{iclass, iidx2(i,1)}{iidx2(i,2)}(1,mod(edgeEnd+1,2)+1);
                                end

                                auxPars = compParValues(edgeEnd, auxVals);

                                if (edgeEnd == edgeIni && parStart > parLast)
                                    iidx2 = iidx2(auxPars > min(parLast, parStart) & auxPars < max(parStart, parLast),:);
                                    auxPars = auxPars(auxPars > min(parLast, parStart) & auxPars < max(parStart, parLast));
                                else
                                    iidx2 = iidx2(auxPars > parLast,:);
                                    auxPars = auxPars(auxPars > parLast);
                                end

                                % Here, jclass is the index in iidx2, not
                                % actually jclass.
                                [~, jclass] = min(auxPars);
                                % If jclass is not empty, this means that
                                % there are fault line components starting
                                % on the current domain boundary. If so, we
                                % take the nearest one. If not, we proceed
                                % to the next domain boundary.
                                if (any(jclass))
                                     jcomp = iidx2(jclass,2);
                                     jclass = iidx2(jclass,1);
                                else
                                    jclass = 0;
                                end
                            end

                            % If there are no other boundary parts leaving
                            % the domain at the current domain boundary, go
                            % to the next corner point, if appropriate.
                            if (jclass <= 0)
                                if (edgeIni == edgeEnd && parStart > parLast)
                                    bfinishedPart = true;
                                else                            
                                    Subdomains{ClassVals(iclass)}{icomp} = ...
                                        [Subdomains{ClassVals(iclass)}{icomp}; ...
                                         CornersNext(edgeEnd,:)];
                                    % No further boundary components on
                                    % current domain boundary: proceed to
                                    % next edge.
                                    edgeEnd = mod(edgeEnd,4) +1;
                                end
                            end

                        % The current component ends inside the domain.
                        else

                            % Search all other boundary parts which start
                            % inside.
                            iidx = (LeftDomainStartAux(iclass, :, :) == 0 & NumPointsSurf(iclass,:, :) > 0);
                            iidx = squeeze(iidx);

                            if (size(iidx,1) == 1)
                                iidx = iidx';
                            end

                            if any(iidx(:))
                                % get the indices of the corresponding boundary
                                % parts
                                [i,j] = ind2sub(size(iidx), find(iidx ==1));
                                iidx2 = zeros(size(i,1), 2);
                                iidx2(:,1) = i;
                                iidx2(:,2) = j;
                                auxVals = zeros(size(iidx2, 1),2);                            

                                % collect all starting points and compute their
                                % distance to the end point of the current part
                                for i = 1: size(iidx2, 1)
                                    auxVals(i,:) = PointSetsSurfaceAux{iclass, iidx2(i,1)}{iidx2(i,2)}(1,:);
                                end
                                auxVals = auxVals - Subdomains{ClassVals(iclass)}{icomp}(end,:);
                                auxVals = auxVals.*auxVals;
                                auxVals(:,1) = auxVals(:,1) + auxVals(:,2);

                                [minDist, jclass] = min(auxVals(:,1));

                                minDist = sqrt(minDist);
                                % if there are candidates which are reasonably
                                % close, continue with them.
                                if (minDist < FaultApproxParams.maxDistForSurfacePoints)
                                % up to here, jclass is the index in iidx2, not
                                % actually jclass 
                                    jcomp = iidx2(jclass,2);
                                    jclass = iidx2(jclass,1);
                                else
                                    bfinishedPart = true;
                                end                                 

                            else
                                if ~any(NumPointsSurf(iclass, :) )
                                    bfinishedSub = true;
                                    bfinishedPart = true;
                                else
                                   warning('Contradicting information about subdomains')
                                   succeeded = false;
                                   return
                                end
                            end
                        end
                    end % There are boundary parts belonging to current subdomain.
                end % while-loop

                %There are no boundary components left: close the polygon.
                switch LeftDomainStartAuxInit

                    % we start at the left domain boundary
                    case(4)

                        switch(edgeEnd)

                            % we end at the bottom
                            case(1)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)  ];

                            % we end at the right:
                            case(2)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)  ];

                            % we end at the top
                            case(3)
                               Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                               ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                               ProblemDescr.Xmin(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)  ];

                            % we end at the left:
                            case(4)
                                numPoints = size(Subdomains{ClassVals(iclass)}{icomp},1);
                                if (Subdomains{ClassVals(iclass)}{icomp}(1,2) > Subdomains{ClassVals(iclass)}{icomp}(numPoints,2))
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmin(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)   ];
                                else
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                        ProblemDescr.Xmin(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)];
                                end
                        end                
                    % we start at the right domain boundary
                    case(2)
                        switch(edgeEnd)
                            % we end at the bottom
                            case(1)
                               Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                               ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                               ProblemDescr.Xmax(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)  ];

                            % we end at the right
                            case(2)
                                numPoints = size(Subdomains{ClassVals(iclass)}{icomp},1);
                                if (Subdomains{ClassVals(iclass)}{icomp}(1,2) < Subdomains{ClassVals(iclass)}{icomp}(numPoints,2))
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmax(1) Subdomains{ClassVals(iclass)}{icomp}(1,2) ];
                                else
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    ProblemDescr.Xmax(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)  ];
                                end

                            % we end at the top
                            case(3)
                               Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                               ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                               ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                               ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                               ProblemDescr.Xmax(1) Subdomains{ClassVals(iclass)}{icomp}(1,2)  ];

                           % we end at the left
                            case(4)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                ProblemDescr.Xmax(1) Subdomains{ClassVals(iclass)}{icomp}(1,2) ];
                        end
                    % we start at the bottom domain boundary
                    case(1)
                        switch(edgeEnd)

                            % we end at the bottom
                            case(1)
                                numPoints = size(Subdomains{ClassVals(iclass)}{icomp},1);
                                if (Subdomains{ClassVals(iclass)}{icomp}(1,1) < Subdomains{ClassVals(iclass)}{icomp}(numPoints,1))
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                    Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmin(2)  ];
                                else
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmin(2) ];
                                end

                            % we end at the right:
                            case(2)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmin(2)];

                            % we end at the top
                            case(3)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmin(2)];

                            % we end at the left:
                            case(4)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmin(2)];
                        end
                    % we start at the top domain boundary
                    case(3)
                        switch(edgeEnd)

                            % we end at the bottom
                            case(1)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmax(2)];

                            % we end at the right:
                            case(2)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmax(2)];

                            % we end at the top
                            case(3)
                                numPoints = size(Subdomains{ClassVals(iclass)}{icomp},1);
                                if (Subdomains{ClassVals(iclass)}{icomp}(1,1) > Subdomains{ClassVals(iclass)}{icomp}(numPoints,1))
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); ...
                                    ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                    ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                    Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmax(2) ];
                                else
                                    Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                    Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmax(2)];
                                end

                            % we end at the left:
                            case(4)
                                Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; ...
                                ProblemDescr.Xmin(1) ProblemDescr.Xmin(2); ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                                ProblemDescr.Xmax(1) ProblemDescr.Xmax(2); ...
                                Subdomains{ClassVals(iclass)}{icomp}(1,1) ProblemDescr.Xmax(2)];
                        end

                    % we start somewhere in the middle 
                    case(0)
                        Subdomains{ClassVals(iclass)}{icomp} = [Subdomains{ClassVals(iclass)}{icomp}; Subdomains{ClassVals(iclass)}{icomp}(1,:)];
                end
                if (bfinishedSub)
                    nparts = icomp;
                    Subdomains{ClassVals(iclass)} = Subdomains{ClassVals(iclass)}(1:nparts);

                    break
                end

            end

            if (any(NumPointsSurf(iclass,:) > 0))
                warning(['There are leftover boundary parts from subdomain ', int2str(ClassVals(iclass)), ', reconstruction of subdomains failed.'])
                succeeded = false;
                return
            end
        end

        % last part: Is class iclass inside or outside the polygon? We
        % know that the fault line components a polygon consists of are
        % ordered such that class iclass is left of the component.
        % If not inside, flip the subdomain aka taking its complement. Of
        % course, this is correct only if the current subdomain is the only
        % one.
        % However, testing for the need of flipping is necessarily heuristic.
        % We choose the point to test with as follow: For two consecutive
        % boundary points, we take the mean and add to that a scaled (intended)
        % inner normal vector. If this point is in fact inside the domain,
        % class iclass is inside.
        % However, testing like this with one point near the boundary and
        % testing with it, it may happen that the point near the boundary
        % which is intended to be inside the domain is in fact outside:
        %              testing point  O
        %                             /\ scaled inner normal vector
        %                             |
        %                           __|____________---------------x
        %  ___________--------------  |                            |
        % x---------------------------o---------------------------x
        % This failure mode is quite rare and does occur merely in the
        % presence of cusps. Making the normal vector too small to prevent this
        % failure mode may lead to numerical instabilities.
        % To minimize the chance for failure, we consider
        % two points instead one one. In the case of contradicting results (one
        % indicates to flip, the other not to flip), we consider a third point
        % and decide by majority.
        % Note that we test with respect to the discrete approximations of a
        % class; therefore, no additional function calls are necessary, only
        % some point-in-polygon-tests.
        for iclass = 1 : nclasses
            for ipart = 1: size(Subdomains{ClassVals(iclass)})

                idx1 = 1;
                % Number of points per subdomain
                idx2 = size(Subdomains{ClassVals(iclass)}{ipart},1);
                idx2 = min(idx2-1, ceil(0.5*idx2));
                auxVec = Subdomains{ClassVals(iclass)}{ipart}([idx1 idx2],:) - Subdomains{ClassVals(iclass)}{ipart}([idx1+1 idx2+1],:);

                testPoint = Subdomains{ClassVals(iclass)}{ipart}([idx1 idx2] ,:) -0.5*auxVec - 0.01*[-auxVec(:, 2),  auxVec(:, 1)];
                inside = inpolygon(testPoint(:,1), testPoint(:,2), Subdomains{ClassVals(iclass)}{ipart}(:,1),Subdomains{ClassVals(iclass)}{ipart}(:,2));

                bflip = false;
                % both two points are outside
                if (all(~inside))
                    bflip = true;

                % in case of ambiguity: let a third point decide
                elseif (any(inside) && any(~inside))
                    idx3 = size(Subdomains{ClassVals(iclass)}{ipart},1);

                    % the following heuristic works if the subdomain consists
                    % of more than three points
                    if (idx3 > 3)
                        idx3 = min(idx3-1, ceil(0.75*idx3));
                        idx3Next = idx3 + 1;
                    else
                        idx3 = 3;
                        idx3Next = 1;
                    end            
                        auxVec = Subdomains{ClassVals(iclass)}{ipart}(idx3,:) - Subdomains{ClassVals(iclass)}{ipart}(idx3Next,:);
                        testPoint = Subdomains{ClassVals(iclass)}{ipart}(idx3,:) -0.5*auxVec - 0.05*[-auxVec(:, 2),  auxVec(:, 1)];
                        inside = inpolygon(testPoint(:,1), testPoint(:,2), Subdomains{ClassVals(iclass)}{ipart}(:,1),Subdomains{ClassVals(iclass)}{ipart}(:,2));
                        bflip = ~inside;
                end

                if (bflip)
                    % Note that polygons in MATLAB do not need to duplicate the
                    % starting points as finishing point; therefore visualizing
                    % it with plot does not lead to a closed contour. This is
                    % not a bug.
                    Subdomains{ClassVals(iclass)}{ipart} = [ProblemDescr.Xmin;  ProblemDescr.Xmax(1) ProblemDescr.Xmin(2); ...
                        ProblemDescr.Xmax; ProblemDescr.Xmin(1) ProblemDescr.Xmax(2); NaN NaN; Subdomains{ClassVals(iclass)}{ipart}];
                end
            end
        end

        % here, any subdomain must consist of at least three points. Check
        % this.
        for iclass = 1: nclasses
            ncomps = size(Subdomains{iclass});
            
            for icomp = 1: ncomps
                if (size(Subdomains{iclass}{icomp}, 1) < 3)
                    succeeded = false;
                    warning(['Component ', int2str(icomp), ...
                             ' of subdomain ', ...
                             int2str(ClassVals(iclass)), ...
                             ' consists of two points only.'])
                    return
                end
            end
        end
        
        succeeded = true;

        % another consistency check: test, that subdomains do not overlap
        % As heuristic, we take test points and test if these are included in
        % more than one part of a subdomains
        testPoints = CreateHaltonSet(200, 2);

        % affine transformation from [0,1]^2 to [Xmin, Xmax]^2
        for idim = 1:2
            testPoints(:,idim) = (ProblemDescr.Xmax(idim) - ProblemDescr.Xmin(idim)).*testPoints(:,idim) + ProblemDescr.Xmin(idim);
        end


        inside = zeros(size(testPoints, 1),1);

        for iclass = 1 : nclasses
            for ipart = 1: size(Subdomains{ClassVals(iclass)})
                inside = inside + inpolygon(testPoints(:,1), testPoints(:,2), Subdomains{ClassVals(iclass)}{ipart}(:,1),Subdomains{ClassVals(iclass)}{ipart}(:,2));
            end
        end

        if (any(inside > 1))
            warning('Overlapping parts of subdomains detected, reconstruction of subdomains failed.')
            succeeded = false;
        end
        
    % end dim == 2
    elseif (dim == 3)

        % evil hack
        for iclass = 1: nclasses
            for jclass = iclass+1:nclasses
                % first part: mesh
                Subdomains{iclass}{jclass}{1} = boundary(PointSetsSurface{iclass, jclass}{1}, 0.5);

                % second part: points
                Subdomains{iclass}{jclass}{2} = 0.5*( + PointSetsSurface{iclass, jclass}{1} + PointSetsSurface{jclass, iclass}{1});
            end
        end
        succeeded = true;
    end
        
    function parVal = compParValues(edgeEnd, coordLast)
        
        switch(edgeEnd)
            case 1
                parVal = coordLast;
            case 2
                parVal = 1 + coordLast;
            case 3
                parVal = 3 - coordLast;
            case 4                
                parVal = 4 - coordLast;
            otherwise
                parVal = 0;
        end
    end

    
end