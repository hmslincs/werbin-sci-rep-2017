function [PL] = CellSelectorRect(PointList)
% CellSelectorRect takes in PointList struct and prompts you to select 
% a region within the image. Once a figure appears use the mouse 
% to define a rectangular region of interest. 
% The method returns only the pnts that are within the region.
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

n= numel(PointList);
pnts=[];
for i=1:n
    pnts=vertcat(pnts,PointList{i}.pnts);
end

img = hist3(pnts(:,1:2),'Edges',{[0:max(pnts(:,1))],[0:max(pnts(:,2))]});
h=figure;
imshow(img,[0,max(img(:)*.01)]);
colormap('hot');
rect = getrect;

%convert rect into a polygon
xi = [rect(1),rect(1),rect(1)+rect(3),rect(1)+rect(3),rect(1)];
yi = [rect(2),rect(2)+rect(4),rect(2)+rect(4),rect(2),rect(2)];



for i=1:n
    pnts = PointList{i}.pnts;
    in = inpolygon(pnts(:,1),pnts(:,2),yi,xi);
    PointList{i}.pnts = pnts(in,:);
    PointList{i}.fullData = PointList{i}.fullData(in,:); 
end

close(h);

PL = PointList;

end
