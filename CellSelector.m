function [CI] = CellSelector(clusterInfo)
% CellSelector takes in clusterInfo and prompts you to select the region 
% containing a cell. Once a figure appears use the mouse to define a 
% polygon. The method returns only the clustuers that are within the region.
%
% Written by Jeffrey Werbin, Ph.D
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

pnts = vertcat(clusterInfo.pnts);
img = hist3(pnts(:,1:2),'Edges',{[0:max(pnts(:,1))],[0:max(pnts(:,2))]});
pnts = vertcat(clusterInfo.ptClusterCenter);
imshow(log(img),[]);
[BW,xi,yi]= roipoly;

in = inpolygon(pnts(:,1),pnts(:,2),yi,xi);

CI= clusterInfo(in);

end
