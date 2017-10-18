function [sucess]=RandomResampleTreeBaggerFunction(comp1,comp2,savedName,nTrees)
% randomly samples an equal number of samples from comp1 and comp2 and runs 
% TreeBagger 100 times
% Makes TreeBagger models with nTrees amd collects the variable importance
% metric from each run.
% 
% Saves the resulting calculations into a file.
%
% Inputs:
%        comp1, a n x k arrray
%        comp2, a m x k array where k is the same as comp2
%    savedName, filename to save results
%       nTrees, number of trees to use for each TreeBagger
%
% Outputs:
%       VarImp, a k x k array that records the how many times a feature had a particular importance rank
%   randVarImp, a k x k array same as above execpt the TreeBagger was trained with randomized labels
%        varSc, a 100 x k array with the OOBPermutedVarDeltaError for each run for each feature
%    randVarSc, same as above but with trained with randomized labels
%
% Written by Jeffrey Werbin, Ph.D.
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

sucess=1;

if nargin <4
   nTrees=50
end

nfeat = size(comp1,2)
VarImp = zeros(nfeat,nfeat);
varSc = zeros(100,nfeat);
randVarImp = zeros(nfeat,nfeat);
randVarSc=varSc;

ParaOps =statset;
ParaOps.UseParallel = true;

%Take an equal number of type 0 and type 1

switched = 0;

if(size(comp1,1)>size(comp2,1))
    %switch comp1 and comp2
    switched=1;
    temp = comp1;
    comp1=comp2;
    comp2=temp;
end


for i = 1:100
    n = fix(numel(comp1(:,1))*0.66);
    ind = randperm(numel(comp1(:,1)),n);
    ind2 = randperm(numel(comp2(:,1)),n);
    
    type = vertcat(zeros(size(ind))',ones(size(ind2))');
    measure = vertcat(comp1(ind,:),comp2(ind2,:));
 b = TreeBagger(nTrees,measure,type,'OOBVarImp','On','MinLeaf',300,'Options',ParaOps);
    [t,idx] = sort(b.OOBPermutedVarDeltaError);
    idx2 = sub2ind([nfeat,nfeat],1:nfeat,idx);
    VarImp(idx2)=VarImp(idx2)+1;
varSc(i,:)=b.OOBPermutedVarDeltaError;

    ind = randperm(numel(type));
b = TreeBagger(nTrees,measure,type(ind),'OOBVarImp','On','MinLeaf',300,'Options',ParaOps);
    [t,idx] = sort(b.OOBPermutedVarDeltaError);
    idx2 = sub2ind([nfeat,nfeat],1:nfeat,idx);
    randVarImp(idx2)=randVarImp(idx2)+1;
randVarSc(i,:)=b.OOBPermutedVarDeltaError;
    
end

['nTrees = ',num2str(nTrees)]
mean(varSc,1)
std(varSc,1)

save(savedName)

end
