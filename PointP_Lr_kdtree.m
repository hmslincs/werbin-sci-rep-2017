function [Lr]=PointP_Lr_kdtree(pos,r)
% PointP_Lr_kdtree takes a point process and calculates Besag's L statistic
% Implemented using KDtreeBallQuery.
%
% Inputs: 
%       pos, reference postion list (nx2)
%         r, a vector of radii to consider
%
% Outputs:
%       Lr, the Lr statistic for 
%
% Jeffrey L. Werbin
% Harvard Medical School
%
% Last Update: 9/30/2013
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

%edge = 1-(4/(3*pi))*((r/L)+(r/W))+((11/(3*pi))-1)*(r^2/(L*W));
edge=1;

%The number of total points
%n = size(posA,1);
%m = size(posB,1);

%pos = [posA ; posB]; %concatinating the two lists allows for processing with one dist matrix

%find boundaries
limX = [max(pos(:,1)), min(pos(:,1))];
limY = [max(pos(:,2)), min(pos(:,2))];
W = limX(1)-limX(2);
L = limY(1)-limY(2);

%to analyze only points > max(rA,rB) from the boundries of the rectangle
%removes edge effect problems

Rmax = max(r);
x = pos(:,1)-limX(2);
y = pos(:,2)-limY(2);

%Here we find all the points in pos A that are at least Rmax away from the
%bounderies. To help prevent edge effects.

%ind = find( (x > Rmax) & (x < W -Rmax) & (y > Rmax) & (y < L-Rmax));
ind = (x > Rmax) & (x < W -Rmax) & (y > Rmax) & (y < L-Rmax);
clear x y;

%keeps only points that meet this criteria as query points
posA= pos(ind,:);

%Uncomputed values are returned as NaNs
Lr = NaN(numel(r),1);

%Average density within the total area considered
lambda = numel(pos(:,1))/(W*L);

%c
timeElap = 0;
for i=1:numel(r)
tic
    %for each r it finds the inPnts that are <r(i) from each query point
    [idx, dist] = KDTreeBallQuery(pos,posA,r(i));
    
    %idx is a cell array with indicies of all posB within r(i) of posA
    % we only need to know the number of points in idx for each point
    temp = cellfun(@(x) numel(x)-1,idx);
    temp = sqrt(sum(temp)/(lambda*pi*(numel(temp)-1)));
    
    
    %out(i,1)= sum(tempB)/(edge*pi*rB^2);
    %out(i,2)= sum(tempA)/(edge*pi*rA^2);
    Lr(i)= temp;
t = toc;
display(['The ',num2str(i),'th r took ',num2str(t),' seconds']);
timeElap = timeElap+t;    
end

timeElap = timeElap/(60*60);

display(['Total time is ',num2str(timeElap),' hours']);


end
    
