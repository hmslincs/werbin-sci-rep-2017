function [clusterInfo] = addToClusterInfo(clusterInfo,points,varargin)
% addToClusterInfo computes extra information for each cluster in clusterInfo 
% For each cluster the x-y coordianate of the points in the cluster, the 
% cluster area in um and the convex hull are added to clusterInfo
%
% Inputs:
%    clusterInfo, struct that is the result of meanshift clustering
%         points, the n x 2 array that was used to generate clusterInfo 
%
% Outputs: 
%    clusterInfo, similar to the above struct with extra infomation added
%
%
% Writen by Jeffrey Werbin, Ph.D. 
% Harvard Medical School, 2014
%

% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
% This file is part of u-track.
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>. 

ip = inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('clusterInfo',@isstruct);
ip.addRequired('points',@(x) (isnumeric(x) & size(x,2) >= 2));

ip.addOptional('DoComposition',true,@islogical);
ip.addOptional('pixelSize',62.81,@isnumeric);

ip.parse(clusterInfo, points,varargin{:});

DoComposition = ip.Results.DoComposition;
pixelSize = ip.Results.pixelSize;

if DoComposition
    n = max(points(:,3));
end


for i=1:numel(clusterInfo)

    pnts = points(clusterInfo(i).ptIdData,:);

    if numel(unique(pnts(:,1))) > 2
    %finds convex hull and area
    [hull,area]= convhull(pnts(:,1:2));
    else
        hull =[];
        area = NaN;
    end
    
    % n is the number of pnts that are merged. Here we store the
    % number of points in this merged cluster from each element of the
    % shift array
    
    if DoComposition
        composition = hist(pnts(:,3),1:n);
        clusterInfo(i).composition = composition;
    end

    
    %Add values to clusterInfo
    clusterInfo(i).pnts = pnts;
    clusterInfo(i).area = area*pixelSize^2;
    clusterInfo(i).hull = hull;
    
end



end
