function [ fullpath ] = wyssFileRipley( file, varargin )
%exchangePaintAnalysis
%   Takes a base directory and creates a list of .mat files that are
%   associate with a single exchangePaint image, then aligns them and
%   preforms some data analysis
%   
%
% Inputs:
%           file, points to a .mat file that contains a PointList of aligned localization movies
%                
%
%Optional:
%
%
%Output: 
%        fullpath: the path to a .mat files containing the results
%
%
%Written by Jeffrey Werbin 
%Harvard Medical School 2014/04/16
%

ip = inputParser;

ip.addRequired('file',@ischar);

ip.addOptional('display',0,@isnumeric);
ip.addOptional('pixelSize',160.5,@isnumeric);


ip.parse(file,varargin{:});

display=ip.Results.display;
pixelSize = ip.Results.pixelSize;

    

%name = [dir(f(end-1)+1:f(end)-1),'_',dir(f(end)+1:end)];
name = [file(1:end-4),'_Ripley.mat'];

load(file,'PointList')

%Used in several places
n = numel(PointList);

%Make continous density estmates using pairwise L function

    ContinuousAnalysis = cell(n);
    Lr = cell([n,1]);
    
    %This loop calculates all the pairwise L(r) functions then modifies it to
    %the H(r) function H(r) = L(r) - r
    %You have to do all as they are not symetric
    
    %these are the radii to use for L(r) cross statistics
    r = 0.2:0.2:5;


    
    for i=1:n
        %calculate Lr (self not cross)
	    Lr{i} = PointP_Lr_kdtree(PointList{i}.pnts,r)-r';
        for j =1:n
           %Computes Besag's L and then renormalizes it to the H statistic
           tmp = PointP_Lr_cross_kdtree(PointList{i}.pnts,PointList{j}.pnts,r);
            
            % This is only a crude way to find the right normalization
            %Const = mean(tmp2(15:end));
            ContinuousAnalysis{i,j} = (tmp)-r';
            
        end
    end


%Appends the additional analysis to fullpath file
    save(name,'ContinuousAnalysis','Lr','r');    
      
    
end


