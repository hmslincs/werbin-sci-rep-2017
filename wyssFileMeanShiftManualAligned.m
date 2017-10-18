function success = wyssFileMeanShiftManualAligned(PointList, name,varargin)
% wyssFileMeanShiftManualAligned takes a cell array containing localization 
% data, merges them into a single matrix, applies meanshift clustering to 
% the data and then saves it to a file with the name +'_MeanShiftClustering.mat'
%
% Inputs:
%    PointList, a cell array with all localization information for an image
%         name, base filename for output
%
% Optional:
%        bandW, the bandwidth in pixels to use for MeanShiftClustering
%
% Outputs:
%        a file with the name name+'_MeanShiftClustering.mat' that contains
%        the clusterInfo, clusterMap and an array containing all 
%        localization positions
%
% Written by Jeffrey Werbin, Ph.D.
% Harvard Medical School 2014
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

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('PointList',@iscell);
ip.addRequired('name',@ischar);

ip.addOptional('bandW',0.3,@isscalar);

ip.parse(PointList,name,varargin{:});

bandW = ip.Results.bandW;


%approx parameters for wyss microscope
Imsize = [256,256];
difLim = 2.1075;

success = 1;

TotalPnts =[];

for i=1:numel(PointList);
   TotalPnts = vertcat(TotalPnts,[PointList{i}.pnts,i*ones(size(PointList{i}.pnts(:,1)))]); 
end



% since pixelSize of wyss system is ~160.5 nm this is ~48 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),bandW,'kernel','flat','flagDebug',true);

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',160.5);

save([name,'_MeanShiftCluster.mat'],'clusterInfo','clusterMap','TotalPnts');

end
