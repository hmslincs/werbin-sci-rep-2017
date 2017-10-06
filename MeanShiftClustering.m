function [clusterInfo, pointToClusterMap, pointTraj] = MeanShiftClustering( ptData, bandwidth, varargin )
% A Fast and Scalable implementation of the mean-shift clustering algorithm using KD-Trees
% 
%     Required Input Arguments:
% 
%                     ptData: specifies the input pointset to be clustered
%                             must be a numPoints x numDimension matrix 
% 
%                  bandwidth: size of the bandwidth to be used in the mean-shift kernels
%                             It must be a scalar value.
%                             
%     Optional Input Arguments:
%         
%                     kernel: specifies the type of kernel used in the non-parametric 
%                             (parzen-window) density function                 
%                              Default: 'gaussian'
%                              Options:                  
%                                 'gaussian' -- uses a gaussian kernel
%                                 'flat' - uses a flat box averaging kernel
%                     
%                     method: specifies the type of mean-shift algorithm to use
%                             Default: 'optimized'
%                             Options:
%                                  'standard': mode-seeking is done for each and every point separately
%                                 'optimized': all nearby points encountered in the mode-seeking path of a point
%                                              are recorded as visited and mode-seeking is not done for them 
%                                              separately. This gives rise to a significant amount of speedup.
%                                  
%                                  
%              maxIterations: specifies the maximum number of iterations/mean-shifts 
%                             that are allowed searching for the mode of a point                        
%                             Default: 500
%         
%         minClusterDistance: specifies the minimum distance between two distinct clusters.
%                             If the distance between two modes is less than this value 
%                             then the two corresponding clusters will be merged.
%                             Default: 2 * bandwidth
%
%              kernelSupport: specifies the maximum distance of points to
%                             include when evaluating the mean-shift kernel.
%                             Default: =2*bandwidth for 'gaussian' and user-specified kernel
%                                      =banwidth for 'flat' kernel;
%                             
%              flagUseKDTree: true/false
%                             specifies whether or not to use the KD-Tree
%                             Default: true (recommended)
%                            
%                  flagDebug: true/false
%                             specifies whether or not to run in debug mode.
%                             In debug mode, a bunch of stuff is/will-be printed 
%                             and plotted for debugging purposes.
%                            
%     Output Arguments:
%                      
%                clusterInfo: A structure array with one entry per cluster found
%                             The structure contains the following elements:
%                                 numPoints - number of points in the cluster
%                                 ptIdData - Ids (row indices of ptData) of 
%                                            the points belonging to the cluster
%                                 ptClusterCenter - coordinates of the cluster center
%                                 
%                                 
%          pointToClusterMap: vector of size equal to the number of points 
%                             in the input point set.
%                             pointToClusterMap(i) - cluster label of point i.   
%                
%         
% Author: Deepak Roy Chittajallu
% 
% References:
% 
% 1) Comaniciu, D. and P. Meer (2002). 
%    "Mean shift: a robust approach toward feature space analysis.", 
%    IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%    24(5): 603-619.
% 

    p = inputParser;
    p.addRequired( 'ptData', @(x) (isnumeric(x) && ndims(x) == 2) );
    p.addRequired( 'bandwidth', @(x) (isscalar(x)) );
    p.addParamValue( 'kernel', 'gaussian', @(x) ( (ischar(x) && ismember(x, {'gaussian', 'flat'})) || (isa(x,'function_handle') && nargin(x) >= 3 && abs(nargout(x)) >= 1)) );
    p.addParamValue( 'method', 'optimized', @(x) ( (ischar(x) && ismember(x, {'standard', 'optimized'})) ) );
    p.addParamValue( 'maxIterations', 500, @(x) (isscalar(x)) );
    p.addParamValue('kernelSupport',[],@(x)(isscalar(x) && x > 0));
    p.addParamValue( 'minClusterDistance', 2 * bandwidth, @(x) (isscalar(x)) );
    p.addParamValue( 'flagUseKDTree', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebug', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( ptData, bandwidth, varargin{:} );
    
    kernelfunc = p.Results.kernel;
    meanShiftMethod = p.Results.method;
    maxIterations = p.Results.maxIterations;
    minClusterDistance = p.Results.minClusterDistance;      
    flagDebug =  p.Results.flagDebug;
    flagUseKDTree = p.Results.flagUseKDTree;
    
    if ischar(p.Results.kernel)
       
        switch p.Results.kernel        
            case 'gaussian'                
                kernelfunc = @update_mean_gaussian_kernel;
            case 'flat'                             
                kernelfunc = @update_mean_flat_kernel;                
        end
                
    else
          kernelfunc = p.Results.kernel;          
    end
    
    if isempty(p.Results.kernelSupport)
        if strcmp(p.Results.kernel,'flat')           
            kernelSupport = bandwidth;
        else
            kernelSupport = 2*bandwidth;
        end
    else
        kernelSupport = p.Results.kernelSupport;
    end

    switch meanShiftMethod
        
        case 'standard'
            
            [clusterInfo, pointToClusterMap,pointTraj] = StandardMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport);
            
        case 'optimized'
            
            [clusterInfo, pointToClusterMap] = OptimizedMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport);            
            pointTraj = {};
    end
    
end

function [clusterInfo, pointToClusterMap] = OptimizedMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport)

    threshConvergence = 1e-3 * bandwidth;    
    [numDataPoints, numDataDims] = size( ptData );           
    
    %% Run mean-shift for each data point and find its mode
    fprintf( 1, '\nMean-shift clustering on a dataset of %d points ...\n', numDataPoints );
    pointClusterVotes = sparse( numDataPoints, 1 );
    flagPointVisited = false( numDataPoints, 1 );
    numPointsProcessed = 0;
    ptIdRandomized = randperm(numDataPoints);
    clusterInfo = [];
    
        % build kd-tree over all points if requested  -- O(n log(n) log(n)) 
        if flagUseKDTree
            
            if flagDebug
                fprintf( 1, '\n\tBuilding KD-Tree ...\n' );        
            end        
            
            kdtree_points = KDTreeSearcher( ptData );  
            
        end
        
        % for each data point climb the non-parametric density function to find its mode
        meanIterationsElapsed = 0;
        if flagDebug
            fprintf( 1, '\n\tFinding the mode for each data point ...\n' );     
        end       
        
        for i = 1:numDataPoints

            curPtId = ptIdRandomized(i);
            
            if flagPointVisited(curPtId)
                continue;
            end
            
            numIterationsElapsed = 0;        
            ptOldMean = ptData(curPtId,:);

            flagPointOnPath = false( numDataPoints, 1 );
            flagPointOnPath( curPtId ) = true;
            
            while true            

                % get nearest points
                if flagUseKDTree
                    [ptIdNearest] = rangesearch(kdtree_points, ptOldMean, kernelSupport);
                    ptIdNearest = ptIdNearest{1};
                else                    
                    ptIdNearest = exhaustive_ball_query( ptData, ptOldMean, kernelSupport);                       
                end
                
                ptNearest = ptData( ptIdNearest, : );                            
                flagPointOnPath( ptIdNearest ) = true;                                
                
                % call kernel to shift mean
                [ptNewMean] = kernelfunc( ptOldMean, ptNearest, bandwidth ); 
                numIterationsElapsed = numIterationsElapsed + 1;

                % check for convergence
                if numIterationsElapsed >= maxIterations || norm(ptNewMean - ptOldMean) < threshConvergence                
                    meanIterationsElapsed = meanIterationsElapsed + numIterationsElapsed;
                    break;
                else
                    ptOldMean = ptNewMean;
                end            

            end               
            
            ptClusterCenter = ptNewMean;
            
            % make a note of all the points that were visited on the path to this mode 
            flagPointVisited = flagPointVisited | flagPointOnPath;            
            
            % check if a close cluster is present already
            blnCloseClusterFound = false;
            closestClusterDist = [];
            closestClusterId = [];
            for cid = 1:numel( clusterInfo )
            
                curClusterDist = norm( ptClusterCenter - clusterInfo(cid).ptClusterCenter );
                if curClusterDist < minClusterDistance 
                    
                    if ~blnCloseClusterFound || (blnCloseClusterFound && curClusterDist < closestClusterDist)                      
                        blnCloseClusterFound = true;
                        closestClusterId = cid;
                        closestClusterDist = curClusterDist;
                    end
                    
                end
                
            end
               
            if ~blnCloseClusterFound
                clusterInfo(end+1).ptClusterCenter = ptClusterCenter;
                pointClusterVotes( flagPointOnPath, numel(clusterInfo) ) = 1;                
            else
                clusterInfo(closestClusterId).ptClusterCenter = 0.5 * (ptClusterCenter + clusterInfo(closestClusterId).ptClusterCenter);
                pointClusterVotes( flagPointOnPath, closestClusterId ) = pointClusterVotes( flagPointOnPath, closestClusterId ) + 1;
            end            
            
            numPointsProcessed = numPointsProcessed + 1;
            
        end        
    
        meanIterationsElapsed = meanIterationsElapsed / numPointsProcessed;
        
        fprintf( 1, '\n\t%d clusters were found ...\n',  numel( clusterInfo ) ); 
        fprintf( 1, '\n\tMode seeking for each data took an average of %d iterations ...\n', round( meanIterationsElapsed ) ); 

        % assign each point to maximum voting cluster
        [ maxvote, pointToClusterMap] = max( pointClusterVotes, [], 2 );
        
    %% post process clusters
    for i = 1:numel( clusterInfo )
       clusterInfo(i).ptIdData = find( pointToClusterMap == i );
       clusterInfo(i).numPoints = numel( clusterInfo(i).ptIdData );
    end
end

function [clusterInfo, pointToClusterMap, pointTraj] = StandardMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport)

    threshConvergence = 1e-3 * bandwidth;    
    [numDataPoints, numDataDims] = size( ptData );           
    
    pointTraj = cell(numDataPoints,1);
    
    %% Run mean-shift for each data point and find its mode
    fprintf( 1, '\nMean-shift clustering on a dataset of %d points ...\n', numDataPoints );
    pointToClusterMap = sparse( numDataPoints, 1 );
    clusterInfo = [];
    
        % build kd-tree over all points if requested  -- O(n log(n) log(n)) 
        if flagUseKDTree
            
            if flagDebug
                fprintf( 1, '\n\tBuilding KD-Tree ...\n' );        
            end        
            kdtree_points = KDTreeSearcher( ptData );  
            
        end        
        
        % for each data point climb the non-parametric density function to find its mode
        meanIterationsElapsed = 0;
        if flagDebug
            fprintf( 1, '\n\tFinding the mode for each data point ...\n' );        
        end       
        
        for i = 1:numDataPoints

            numIterationsElapsed = 0;        
            ptOldMean = ptData(i,:);
            
            pointTraj{i} = nan(maxIterations,numDataDims);
            pointTraj{i}(1,:) = ptOldMean;

            while true            

                % get nearest points
                if flagUseKDTree
                    [ptIdNearest] = rangesearch(kdtree_points, ptOldMean, kernelSupport);
                    ptIdNearest = ptIdNearest{1};
                else                    
                    ptIdNearest = exhaustive_ball_query( ptData, ptOldMean, kernelSupport );                       
                end
                ptNearest = ptData( ptIdNearest, : );
                
                % call kernel to shift mean
                [ptNewMean] = kernelfunc( ptOldMean, ptNearest, bandwidth ); 
                numIterationsElapsed = numIterationsElapsed + 1;
                
                pointTraj{i}(numIterationsElapsed+1,:) = ptNewMean;

                % check for convergence
                if numIterationsElapsed >= maxIterations || norm(ptNewMean - ptOldMean) < threshConvergence                
                    meanIterationsElapsed = meanIterationsElapsed + numIterationsElapsed;
                    break;
                else
                    ptOldMean = ptNewMean;
                end            

            end               

            ptClusterCenter = ptNewMean;
            pointTraj{i} = pointTraj{i}(1:numIterationsElapsed+1,:);
            
            % get point density around the cluster center
            if flagUseKDTree
                [ptIdNearest] = rangesearch(kdtree_points, ptOldMean, kernelSupport);
                ptIdNearest = ptIdNearest{1};
            else                    
                ptIdNearest = exhaustive_ball_query( ptData, ptClusterCenter, kernelSupport );                       
            end            
            curClusterCenterDensity = numel( ptIdNearest );
            
            % check if cluster is present already
            blnCloseClusterFound = false;
            for cid = 1:numel( clusterInfo )
            
                if norm( ptClusterCenter - clusterInfo(cid).ptClusterCenter ) < minClusterDistance
                    
                    blnCloseClusterFound = true;
                    if curClusterCenterDensity > clusterInfo(cid).clusterCenterDensity
                        clusterInfo(cid).ptClusterCenter = ptClusterCenter;
                        clusterInfo(cid).clusterCenterDensity = curClusterCenterDensity;
                    end
                    pointToClusterMap(i) = cid;
                    
                end
                
            end
               
            if ~blnCloseClusterFound
                clusterInfo(end+1).ptClusterCenter = ptClusterCenter;
                clusterInfo(end).clusterCenterDensity = curClusterCenterDensity;
                pointToClusterMap(i) = numel( clusterInfo );
            end            
            
        end        
    
        meanIterationsElapsed = meanIterationsElapsed / numDataPoints;
        
        fprintf( 1, '\n\t%d clusters were found ...\n',  numel( clusterInfo ) ); 
        fprintf( 1, '\n\tMode seeking for each data took an average of %d iterations ...\n', round( meanIterationsElapsed ) ); 
        
    %% post process clusters
    for i = 1:numel( clusterInfo )
       clusterInfo(i).ptIdData = find( pointToClusterMap == i );
       clusterInfo(i).numPoints = numel( clusterInfo(i).ptIdData );
    end

end


function [ptNewMean] = update_mean_gaussian_kernel( ptOldMean, ptNearest, bandwidth )

    numNeighPoints = size( ptNearest, 1 );
    sqdistances = sum( (ptNearest - repmat(ptOldMean, numNeighPoints, 1) ).^2, 2 );
    weights = exp( -0.5 * sqdistances / bandwidth^2 );    
    ptNewMean = (weights' * ptNearest) / sum(weights);
    
end

function [ptNewMean] = update_mean_flat_kernel( ptOldMean, ptNearest, bandwidth )

    ptNewMean = mean( ptNearest, 1 );
    
end

function [ ptIdNearest ] = exhaustive_ball_query( ptData, ptPoint, bandwidth )

    numDataPoints = size( ptData, 1 );
    sqDistances = sum( ( ptData - repmat( ptPoint, numDataPoints, 1 ) ).^2, 2 );
    ptIdNearest = find( sqDistances < bandwidth^2 );      
    
end
