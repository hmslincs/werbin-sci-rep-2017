function success = wyssFileMeanShiftManualAligned(PointList, name,varargin)
% Takes a cell array containing localization data, merges them into a single
% matrix, applies meanshift clustering to the data and then saves it to 
% a file with the name +'_MeanShiftClustering.mat'
% 


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
