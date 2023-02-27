% This function provides the 2d-specific part of the algorithm. For
% details, we refer to "Detecting and approximating decision boundaries in
% low dimensional spaces" (http://arxiv.org/abs/2302.08179).
%
% Input:
% - PointSetsSurface: (nclasses x nclasses)-structure of structures
%   containing the components of the point sets between subsets.
% - NumPointsSurf (nclasses x nclasses)-structure of arrays containing
%   the number of points in the point sets
% - nclasses: number of different classes
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
%   containing the components of the point sets between subsets.
% - NumPointsSurf (nclasses x nclasses)-structure of arrays containing
%   the number of points in the point sets
% - LeftDomainStart: In 2D, this (nclasses x nclasses x #components)-array
%   codes information on whether a fault line starts on the domain
%   boundary (value > 0, we refer to the code below for ASCII-art
%   explaining how) or inside the domain (value = 0).
% - LeftDomainEnd: The same as LeftDomainStart, but for ending of fault
%   lines
% - bsuccessful: flag, if the fault lines have been processed successfully

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointSetsSurface, NumPointsSurf, LeftDomainStart, ...
          LeftDomainEnd, bsuccessful] = ...
    getTripletsNearFault2D(PointSetsSurface, NumPointsSurf, nclasses, ...
                           ClassVals, FaultApproxParams, ProblemDescr)

    global ncalls
    global ExtendedStats;

    bsuccessful = false;

    % desired maximum distance of a point on the fault line to the next one
    maxDistForSurfacePoints = FaultApproxParams.maxDistForSurfacePoints;
    
    NumCompsPerFaultLine = zeros(nclasses);
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses
            
            if ProblemDescr.verboseMode
                disp(['-- fill gaps on boundary between classes ' ...
                      int2str(ClassVals(iclass)) ' and ' int2str(ClassVals(jclass))])
            end

            % Up to this point, there is only one component per fault
            % line.
            if (NumPointsSurf{iclass, jclass}(1) > 0)

                % Fill possible gaps in the fault line by adding more
                % points in between consecutive points if necessary.
                
                itry = 1;
                
                % We need to repeat this process: It may happen that due
                % to filling some gaps in the first pass, the assumed
                % shape of the fault line changes substantially due to
                % additional infomation at hand.
                % In this case, new gaps may arise which then have to
                % be filled again. We repeat that process until no more
                % gaps are found and filled or the maximal number of
                % filling passes has been reached.
                pointsAdded = true;
                while pointsAdded && itry <= FaultApproxParams.maxTrialsForFillingGaps
                    [PointSetsSurface{iclass, jclass}{1}, ...
                     PointSetsSurface{jclass, iclass}{1}, ...
                     NumPointsSurf{iclass, jclass}(1), ...
                     pointsAdded, bsuccessfulFillGaps] = ...
                        fill2D(PointSetsSurface{iclass, jclass}{1}, ...
                               PointSetsSurface{jclass, iclass}{1}, ...
                               iclass, jclass, ClassVals, ...
                               FaultApproxParams, ProblemDescr, [0,0,1]);

                    if ~bsuccessfulFillGaps
                        warning(['GetPointsNearSurface failed for classes ' ...
                                 int2str(ClassVals(iclass)) ' and ' ...
                                 int2str(ClassVals(jclass)) '.'])

                        % dummy return values
                        LeftDomainStart = 0;
                        LeftDomainEnd = 0;
                        return
                    end
                                        
                    % In fact, points have been added on both sides of the
                    % fault line simultaneously, forming triplets.
                    NumPointsSurf{jclass, iclass}(1) = NumPointsSurf{iclass, jclass}(1);

                    % sort again
                    [IdxPointsSurfOrdered, sortingSuccessful] = ...
                        sortPointsOnFaultLine(PointSetsSurface{iclass, jclass}{1}, ...
                                              1, ProblemDescr, FaultApproxParams);
                    PointSetsSurface{iclass, jclass}{1} = PointSetsSurface{iclass, jclass}{1}(IdxPointsSurfOrdered{1},:);
                    PointSetsSurface{jclass, iclass}{1} = PointSetsSurface{jclass, iclass}{1}(IdxPointsSurfOrdered{1},:);

                    if ~sortingSuccessful
                        warning(['GetPointsNearSurface failed for classes ' ...
                                 int2str(ClassVals(iclass)) ' and ' ...
                                 int2str(ClassVals(jclass)) '.'])

                        % dummy return values
                        LeftDomainStart = 0;
                        LeftDomainEnd = 0;
                        return
                    end
                    
                    itry = itry + 1;
                end
                
                
                % Test, if the fault line between class iclass and
                % jclass consists in fact of several unconnected
                % components.
                [PointSetsSurface{iclass, jclass}, ...
                 NumCompsPerFaultLine(iclass, jclass), ...
                 NumPointsSurf{iclass, jclass}] = ...
                    divideIntoComponents(PointSetsSurface{iclass, jclass}{1}, ...
                                         maxDistForSurfacePoints);
                
                for icomp = 1:NumCompsPerFaultLine(iclass, jclass)
                    noIntersection = ~selfIntersection(PointSetsSurface{iclass,jclass}{icomp});
                    
                    if (~noIntersection)
                        warning(['Self-intersecting boundary components between classes ' ...
                            int2str(ClassVals(iclass)) ' and ' int2str(ClassVals(jclass)) ...
                            ' detected. Skip computation.'])

                        LeftDomainStart = 0;
                        LeftDomainEnd = 0;

                        return
                    end
                end
                
                % Repeat sorting for the counterparts on the other side
                % of the fault line.
                NumCompsPerFaultLine(jclass, iclass) = NumCompsPerFaultLine(iclass, jclass);

                % Hack: as iclass < jclass, NumPointsSurf{jclass, iclass}
                % has not been updated yet and therefore contains the
                % number of points in S_i,j prior to separation in
                % components.
                % Therefore, this entry contains the total number of points
                % in S_i,j.
                numPointsTotal = NumPointsSurf{jclass, iclass}(1);
                
                NumPointsSurf{jclass, iclass} = NumPointsSurf{iclass, jclass};

                PointsTemp = cell(NumCompsPerFaultLine(jclass, iclass),1);
                istart = 1;

                iend = size(PointSetsSurface{iclass, jclass}{1},1);
                for icomp = 1: NumCompsPerFaultLine(jclass, iclass)-1
                    PointsTemp{icomp} = PointSetsSurface{jclass, iclass}{1}(istart:iend,:);
                    istart = iend + 1;
                    iend = istart + size(PointSetsSurface{iclass, jclass}{icomp+1},1) - 1;
                end
                PointsTemp{NumCompsPerFaultLine(jclass, iclass)} = ...
                    PointSetsSurface{jclass, iclass}{1}(istart:numPointsTotal,:);
                PointSetsSurface{jclass, iclass} = PointsTemp;
            end
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_filling_gaps';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end


    % maximal number of components per fault line
    maxNumComps = max(NumCompsPerFaultLine, [], 'all');

    LeftDomainStart = -ones(nclasses, nclasses, maxNumComps);
    LeftDomainEnd = -ones(nclasses, nclasses, maxNumComps);
    
    % Extend all components of a fault line until near to their
    % true start and end points.
    % When extended, perform adaptive refinement/coarsening.
    %
    % indices of the domain edges:
    %        y /\
    %          |       3
    %          |_______________ 
    %          |               |
    %          |               |
    %          |               |
    %        4 |               | 2
    %          |               |
    %          |               |
    %          |_______________|__________ x
    %                  1
    %
    if ProblemDescr.verboseMode
        disp('-- expand boundaries')
    end

    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses
            if (NumPointsSurf{iclass, jclass}(1) > 0)
                for icomp = 1: NumCompsPerFaultLine(iclass, jclass)
                    
                    numPointsComp = NumPointsSurf{iclass, jclass}(icomp);
                    IdxPointsSurfOrdered = 1:numPointsComp;
                    % Try to extrapolate the fault line by fitting a
                    % poynomial curve to the first and the last points in
                    % the point set.

                    % If the component consists of one point only: find
                    % more by scattering.
                    if numPointsComp == 1
                        SinglePoint = PointSetsSurface{iclass, jclass}{icomp};
                        AuxVec = [-1.5, -1; 0.5, -1; 0.5,1; -1.5, 1];
                        
                        pairsFound = false;
                        for itry = 1: 3
                            % scatter four points around the single one
                            Xtest = SinglePoint + 1/2^itry*FaultApproxParams.maxDistForSurfacePoints*AuxVec;
                            classAux = computeClassification(Xtest, ProblemDescr);

                            % points in class iclass
                            Xtest1 = Xtest(classAux == ClassVals(iclass),:);

                            % points in class jclass
                            Xtest2 = Xtest(classAux == ClassVals(jclass),:);

                            % Build all combinations between these points.
                            npointsi = size(Xtest1,1);
                            npointsj = size(Xtest2,1);
                            ncombis = npointsi * npointsj;
                            Pointsiclass = zeros(ncombis,2);
                            Pointsjclass = zeros(ncombis,2);
                            for i = 1: npointsi
                                Pointsjclass(1 + (i-1)*npointsj: i*npointsj,:) = Xtest2;
                            end

                            for i = 1: npointsi
                                Pointsiclass(1 + (i-1)*npointsj: i*npointsj,:) = ...
                                    ones(npointsj,1)*Xtest1(i,:);
                            end

                            for i = 1: ncombis
                                [PointLeft, PointRight, finished] = ...
                                    singleTripletByBisection(Pointsiclass(i,:), ...
                                                             Pointsjclass(i,:), ...
                                                             ClassVals(iclass), ...
                                                             ClassVals(jclass), ...
                                                             ProblemDescr, ...
                                                             FaultApproxParams);

                                if (finished)
                                    NumPointsSurf{iclass, jclass}(icomp) = NumPointsSurf{iclass, jclass}(icomp) + 1;
                                    NumPointsSurf{jclass, iclass}(icomp) = NumPointsSurf{jclass, iclass}(icomp) + 1;

                                    PointSetsSurface{iclass, jclass}{icomp}(NumPointsSurf{iclass, jclass}(icomp),:) = PointLeft;
                                    PointSetsSurface{jclass, iclass}{icomp}(NumPointsSurf{jclass, iclass}(icomp),:) = PointRight;
                                    pairsFound = true;
                                end
                            end
                            if pairsFound
                                break
                            end
                        end
                        numPointsComp = NumPointsSurf{iclass, jclass}(icomp);
                        [IdxPointsSurfOrdered, sortingSuccessful] = ...
                            sortPointsOnFaultLine(PointSetsSurface{iclass, ...
                                                  jclass}{icomp}, 1, ...
                                                  ProblemDescr, ...
                                                  FaultApproxParams);
                        IdxPointsSurfOrdered = IdxPointsSurfOrdered{1};
                    end
                    
                    if numPointsComp == 1
                        warning(['There is only one point on the boundary between classes ', int2str(ClassVals(iclass)), ...
                            ' and ', int2str(ClassVals(jclass)), '. Searching more by scattering failed.'])
                        warning('Finding subdomains terminated.')

                        % dummy return values
                        LeftDomainStart = 0;
                        LeftDomainEnd = 0;
                        return
                    end
                    
                    % beginning of fault line
                    [IdxPointsSurfOrderedNew, LeftDomainStart, NumPointsSurf{iclass,jclass}(icomp), ...
                        PointSetsSurface{iclass,jclass}{icomp}, PointSetsSurface{jclass,iclass}{icomp}, ~, ~] = ...
                        expand2D(LeftDomainStart, IdxPointsSurfOrdered, iclass, jclass,  icomp, ...
                        ClassVals, PointSetsSurface{iclass,jclass}{icomp}, ...
                        PointSetsSurface{jclass,iclass}{icomp}, [0,0,1], ProblemDescr,FaultApproxParams, 0);
                    IdxPointsSurfOrdered = IdxPointsSurfOrderedNew;
                    NumPointsSurf{jclass,iclass}(icomp) = NumPointsSurf{iclass,jclass}(icomp);
                    
                    PointSetsSurface{iclass, jclass}{icomp} = ...
                        PointSetsSurface{iclass, jclass}{icomp}(IdxPointsSurfOrdered, :);
                    PointSetsSurface{jclass, iclass}{icomp} = ...
                        PointSetsSurface{jclass, iclass}{icomp}(IdxPointsSurfOrdered, :);

                    IdxPointsSurfOrdered = 1: NumPointsSurf{iclass, jclass}(icomp);
                    
                    % end of fault line
                    [IdxPointsSurfOrderedNew, LeftDomainEnd, NumPointsSurf{iclass,jclass}(icomp), ...
                        PointSetsSurface{iclass,jclass}{icomp}, PointSetsSurface{jclass,iclass}{icomp}, reSort, ~] = ...
                        expand2D(LeftDomainEnd, IdxPointsSurfOrdered, iclass, jclass, icomp, ...
                        ClassVals, PointSetsSurface{iclass,jclass}{icomp}, ...
                        PointSetsSurface{jclass,iclass}{icomp}, [0,0,1], ProblemDescr, FaultApproxParams, 1);
                    IdxPointsSurfOrdered = IdxPointsSurfOrderedNew;
                    NumPointsSurf{jclass,iclass}(icomp) = NumPointsSurf{iclass,jclass}(icomp);

                    PointSetsSurface{iclass, jclass}{icomp} = ...
                    PointSetsSurface{iclass, jclass}{icomp}(IdxPointsSurfOrdered, :);
                    PointSetsSurface{jclass, iclass}{icomp} = ...
                        PointSetsSurface{jclass, iclass}{icomp}(IdxPointsSurfOrdered, :);
                    
                    % Resort the points on the fault line if indicated.
                    if reSort
                        [IdxPointsSurfOrdered, sortingSuccessful] = ...
                            sortPointsOnFaultLine(PointSetsSurface{iclass, jclass}{icomp}, ...
                                                  1, ProblemDescr, FaultApproxParams);
                        PointSetsSurface{iclass, jclass}{icomp} = ...
                            PointSetsSurface{iclass, jclass}{icomp}(IdxPointsSurfOrdered{1}, :);
                        PointSetsSurface{jclass, iclass}{icomp} = ...
                            PointSetsSurface{jclass, iclass}{icomp}(IdxPointsSurfOrdered{1}, :);
                    
                        if ~sortingSuccessful
                            warning (['resorting of points during expanding the boundary between classes ' ...
                                int2str(ClassVals(iclass)) ' and ' int2str(ClassVals(jclass)) ' failed.'])
                            % dummy return values
                            LeftDomainStart = 0;
                            LeftDomainEnd = 0;
                            return
                        end
                    end
                
                    % Remove duplicates: Duplicates should actually not
                    % occur, but they can due to wrong sorting which
                    % remained unnoticed.
                    [PointSetsSurface{iclass, jclass}{icomp}, PointSetsSurface{jclass, iclass}{icomp}] = ...
                        removeDuplicates(PointSetsSurface{iclass, jclass}{icomp}, ...
                                         PointSetsSurface{jclass, iclass}{icomp}, ...
                                         FaultApproxParams.eps);
                    
                    NumPointsSurf{iclass, jclass}(icomp) = size(PointSetsSurface{iclass, jclass}{icomp}, 1);
                    NumPointsSurf{jclass, iclass}(icomp) = size(PointSetsSurface{jclass, iclass}{icomp}, 1);
                end

                % Consistency check: test, if the different components
                % intersect. It may happen that due to failed sorting,
                % some single point or so was errornously considered a
                % separate boundary component. Now, additional points
                % have been found. If these different boundary components
                % are really different ones, they do not intersect.
                if NumCompsPerFaultLine(iclass,jclass) > 1
                    for icomp = 1: NumCompsPerFaultLine(iclass,jclass)
                        
                        jcomp = icomp+1;
                        % A while-loop is better here, as
                        % NumCompsPerFaultLine may decrease by merging,
                        % which is not reflected in an ordinary for-loop.
                        while jcomp <= NumCompsPerFaultLine(iclass,jclass)
                            doNotIntersect = ...
                                ~polyLinesIntersect(PointSetsSurface{iclass, jclass}{icomp}, ...
                                                    PointSetsSurface{iclass, jclass}{jcomp});
                            pointOnLine = ...
                                isOnPolyLine(PointSetsSurface{iclass, jclass}{icomp}, ...
                                             PointSetsSurface{iclass, jclass}{jcomp});

                            % If two components intersect, they must in
                            % fact be one component. Therefore, we merge
                            % them.
                            if (~doNotIntersect || pointOnLine)
                                
                                % Maybe these two components are in fact
                                % the same: try to merge them.
                                test = [PointSetsSurface{iclass, jclass}{icomp}; ...
                                        PointSetsSurface{iclass, jclass}{jcomp}];
                                [IdxPointsSurfOrdered, sortingSuccessful] = ...
                                    sortPointsOnFaultLine(test, 1, ...
                                                          ProblemDescr, ...
                                                          FaultApproxParams);
                                
                                if (sortingSuccessful)
                                    test = test(IdxPointsSurfOrdered{1},:);
                                else
                                    warning(['Two components of the boundary between classes ' int2str(ClassVals(iclass)) ' and ' ...
                                        int2str(ClassVals(jclass)) ' were found to intersect and their merge failed.']);
                                    return;
                                end
                                
                                % Test, if the new boundary component does
                                % not intersect itself.
                                doNotIntersect = ~selfIntersection(test);
                                if (doNotIntersect)
                                    PointSetsSurface{iclass, jclass}{icomp} = test;
                                    PointSetsSurface{iclass, jclass}(jcomp) = [];

                                    NumPointsSurf{iclass, jclass}(icomp) = ...
                                        NumPointsSurf{iclass, jclass}(icomp) + NumPointsSurf{iclass, jclass}(jcomp);
                                    NumPointsSurf{iclass, jclass}(jcomp) = [];

                                    NumPointsSurf{jclass, iclass}(icomp) = ...
                                        NumPointsSurf{jclass, iclass}(icomp) + NumPointsSurf{jclass, iclass}(jcomp);
                                    NumPointsSurf{jclass, iclass}(jcomp) = [];

                                    PointSetsSurface{jclass, iclass}{icomp} = ...
                                        [PointSetsSurface{jclass, iclass}{icomp}; PointSetsSurface{jclass, iclass}{jcomp}];
                                    PointSetsSurface{jclass, iclass}(jcomp) = [];
                                    
                                    NumCompsPerFaultLine(iclass,jclass) = NumCompsPerFaultLine(iclass,jclass) -1;
                                    NumCompsPerFaultLine(jclass,iclass) = NumCompsPerFaultLine(jclass,iclass) -1;
                                    
                                    % Reset LeftDomainStart and
                                    % LeftDomainEnd, as the corresponding
                                    % component is gone.
                                    LeftDomainStart(iclass, jclass, jcomp:end-1) = ...
                                        LeftDomainStart(iclass, jclass, jcomp+1:end);
                                    LeftDomainEnd(iclass, jclass, jcomp:end-1) = ...
                                        LeftDomainEnd(iclass, jclass, jcomp+1:end);
                                    LeftDomainEnd(iclass, jclass, end) = -1;
                                    LeftDomainStart(iclass, jclass, jcomp) = -1;
                                    
                                    % maximal number of components per
                                    % fault line
                                    maxNumCompsNew = max(NumCompsPerFaultLine, [], 'all');
                                    
                                    % If the maximal number of components
                                    % per fault line decreased due to
                                    % merging, shorten the arrays
                                    % accordingly.
                                    if (maxNumCompsNew < maxNumComps)
                                        LeftDomainStart = LeftDomainStart(:,:,1:maxNumCompsNew);
                                        LeftDomainEnd = LeftDomainEnd(:,:,1:maxNumCompsNew);
                                        maxNumComps = maxNumCompsNew;
                                    end
                                    
                                else
                                    warning(['Two components of the boundary between classes ' ...
                                             int2str(ClassVals(iclass))...
                                            ' and ' int2str(ClassVals(jclass)) ...
                                            ' were found to intersect. However, their merge']);
                                    warning('led to another intersecting boundary components. Stop computation.')
                                    return;
                                end
                                
                            else
                                jcomp = jcomp+1;
                            end
                        end
                    end
                end
            end
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_prol_lines';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end

    % adaptive refinement and coarsening according to estimated curvature
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses

            if ProblemDescr.verboseMode
                disp(['-- adaptive refinement on boundary between classes ' ...
                    int2str(ClassVals(iclass)) ' and ' int2str(ClassVals(jclass))])
            end

            if (NumPointsSurf{iclass, jclass}(1) > 0)
                for icomp = 1:NumCompsPerFaultLine(iclass,jclass)
                   [PointSetsSurface, NumPointsSurf] = adapt2D(PointSetsSurface, NumPointsSurf, ...
                       iclass, jclass, icomp, ProblemDescr, FaultApproxParams, ClassVals);
                end
            end
        end
    end

    if ProblemDescr.extendedStats
        ExtendedStats.pos_in_code{end+1} = 'after_adaptive_ref';
        ExtendedStats.ncalls{end+1} = ncalls(2);
        ExtendedStats.PointSetsSurf{end+1} = PointSetsSurface;
        ExtendedStats.nPointsSurf{end+1} = NumPointsSurf;
    end
    
    % If necessary, reverse the order of fault line components.
    for iclass = 1:nclasses
        for jclass = iclass+1:nclasses
            for icomp = 1: NumCompsPerFaultLine(iclass, jclass)
                if (NumPointsSurf{iclass, jclass}(icomp) > 1)
    
                    % Sort the points near the fault line such that the
                    % subdomain iclass is right to the line.

                    % Take a point from somewhere in the middle of the
                    % line.
                    iidx = max(2, round(NumPointsSurf{iclass, jclass}(icomp)/2));

                    % vector pointing from point with index iidx to
                    % iidx-1
                    AuxVec = PointSetsSurface{iclass, jclass}{icomp}(iidx,:) - ...
                             PointSetsSurface{iclass, jclass}{icomp}(iidx-1,:);
                    NormalVec = FaultApproxParams.alpha*[AuxVec(2), -AuxVec(1)];
                    
                    % The points are known up to 2*abstolBisection
                    % only. Therefore, a normal vector smaller than
                    % this does not make sense.
                    % If the boundary points are very close, our simple
                    % estimation of the normal vector may be
                    % inaccurate, such that we should elongate it even
                    % a little more to be safe that in case of wrong
                    % ordering, the auxiliary point is really beyond
                    % the boundary of the subdomain. This leads to the
                    % factor of 3.
                    if (norm(NormalVec) < 3*FaultApproxParams.abstolBisection)
                        NormalVec = NormalVec/norm(NormalVec)*3*FaultApproxParams.abstolBisection;

                    % The normal vector can be overly large (this
                    % can happen e.g. after adaptive refinement and
                    % coarsening of a straight line, which is represented
                    % by a few points only). Then, we shorten it to a size
                    % of maxDistForSurfacePoints.
                    elseif (norm(NormalVec) > maxDistForSurfacePoints)
                        NormalVec = NormalVec/norm(NormalVec)*maxDistForSurfacePoints;
                    end

                    % The aux point should be in class ClassVals(iclass).
                    % If not, either the point is badly chosen or the
                    % order is wrong. 
                    % Special treatment is required if the point with
                    % index iidx is the end of the boundary component.
                    % This may happen if e.g. the boundary component
                    % consists of two points only.
                    % Above, we excluded that iidx is the first point
                    % on the boundary component.
                    % The first or the last point on a boundary component
                    % is unsuitable for our purpose.

                    % Point is somewhere in the middle: just take that
                    % point. Note that iidx >=2, such that
                    % NumPointsSurf{iclass, jclass}(icomp) >= 3
                    if (iidx < NumPointsSurf{iclass, jclass}(icomp))
                        auxPoint = PointSetsSurface{iclass, jclass}{icomp}(iidx,:) - NormalVec;

                    % Very last point, more than two points on the
                    % component: take the predecessor.
                    elseif (iidx == NumPointsSurf{iclass, jclass}(icomp) && ...
                            NumPointsSurf{iclass, jclass}(icomp) > 2)
                        auxPoint = PointSetsSurface{iclass, jclass}{icomp}(iidx-1,:) - NormalVec;

                    % Very last point, and just two points on the
                    % component: compute an auxiliary point based on
                    % the mean of the two points.
                    else
                        auxPoint = PointSetsSurface{iclass, jclass}{icomp}(iidx-1,:) + 0.5*AuxVec - NormalVec;
                    end
                    classAux = computeClassification(auxPoint, ProblemDescr);
                    if (classAux == ClassVals(jclass))
                        PointSetsSurface{iclass, jclass}{icomp} = flip(PointSetsSurface{iclass, jclass}{icomp});
                        
                        % Start point becomes end point.
                        aux = LeftDomainEnd(iclass, jclass, icomp);
                        LeftDomainEnd(iclass, jclass, icomp) = LeftDomainStart(iclass, jclass, icomp);
                        LeftDomainStart(iclass, jclass, icomp) = aux;
                    elseif(classAux == ClassVals(iclass))
                        PointSetsSurface{jclass, iclass}{icomp} = flip(PointSetsSurface{jclass, iclass}{icomp});

                    % It might happen that the auxiliary point is
                    % outside the domain or inside a entirely different
                    % subdomain. In this case, we reverse it
                    % and try again. If it still fails, we give up.
                    else
                        % Point is somewhere in the middle: just take that
                        % point. Note that iidx >=2, such that
                        % NumPointsSurf{iclass, jclass}(icomp) >= 3.
                        if (iidx < NumPointsSurf{iclass, jclass}(icomp))
                            auxPoint = PointSetsSurface{iclass, jclass}{icomp}(iidx,:) + NormalVec;
                        % Very last point, more than two points on the
                        % component: take the predecessor.
                        elseif (iidx == NumPointsSurf{iclass, jclass}(icomp) && ...
                                NumPointsSurf{iclass, jclass}(icomp) > 2)
                            auxPoint = PointSetsSurface{iclass, jclass}{icomp}(iidx-1,:) + NormalVec;
                        % Very last point, and just two points on the
                        % component: compute an auxiliary point based on
                        % the mean of the two points.
                        else
                            auxPoint = PointSetsSurface{iclass, jclass}{icomp}(iidx-1,:) + 0.5*AuxVec + NormalVec;
                        end

                        classAux2 = computeClassification(auxPoint, ProblemDescr);                        

                        if (classAux2 == ClassVals(iclass))
                            PointSetsSurface{iclass, jclass}{icomp} = flip(PointSetsSurface{iclass, jclass}{icomp});

                            % Start point becomes end point.
                            aux = LeftDomainEnd(iclass, jclass, icomp);
                            LeftDomainEnd(iclass, jclass, icomp) = LeftDomainStart(iclass, jclass, icomp);
                            LeftDomainStart(iclass, jclass, icomp) = aux;
                        elseif(classAux2 == ClassVals(jclass))
                            PointSetsSurface{jclass, iclass}{icomp} = flip(PointSetsSurface{jclass, iclass}{icomp});
                        else
                            warning(['Determination of orientation failed for the boundary between classes ' ...
                                num2str(ClassVals(iclass)) ' and ' num2str(ClassVals(jclass)) '.']);
                            return;
                        end
                    end
                end
            end
        end
    end
    bsuccessful = true;
end