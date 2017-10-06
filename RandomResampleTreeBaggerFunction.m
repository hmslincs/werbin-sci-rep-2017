function [sucess]=RandomResampleTreeBaggerFunction(comp1,comp2,savedName,nTrees)
% randomly samples type = 0 and type = 1 and runs TreeBagger 50 times
% uses 50 trees
% returns of 0.05 confidence levels by permutation
%
% 
%

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
