# Bottleneck distance 

You can download jplex toolbox [here](https://www.math.colostate.edu/~adams/jplex/files/PlexMatlabTutorial.pdf) for estimating holes in a simplicial complx. 

![toyexample_topologicalspace](https://user-images.githubusercontent.com/54297018/63508330-b8165580-c514-11e9-98b6-570b28aa84f3.png)

We randomly select 100 points in each topological space, and save data points in 'x.mat'. 

```Matlab 
% load data set 
load 'x.mat'; 

% Number of nodes 
p = size(x,1);
% Index of the upper triangular part of a connectivity matrix 
ind_triu = find(triu(ones(p,p),1));
% Number of edges 
q = length(ind_triu);

% Estimate Euclidean distance between point clouds. 
d = [];
for g = 1:4, % number of topological spaces 
    for i = 1:p,
        for j = i+1:p,
            d(i,j,g) = sqrt(sum((x(i,:,g) - x(j,:,g)).^2));
            d(j,i,g) = d(i,j,g);
        end
    end
end

% Do jplex for estimating the birth and death of the 0- and 1-dimensional holes in four sets of point cloud data 
% Parameters for jplex
max_dimension = 2;
max_filtration_value = 2500;
num_divisions = 2500;

barcode0_D = []; barcode1_D = [];
for g = 1:4, 
    D = d(:,:,g);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % This Hyekyoung Lee's own way to estimate holes using jplex
    % You don't have to follow. 
    % Perform.
    [sval,sind] = sort(D(ind_triu),'ascend');
    tmp = zeros(q,1);
    tmp(sind) = [1:q];
    distances = zeros(p,p);
    distances(ind_triu) = tmp;
    distances = distances + distances';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do javaplex.
    m_space = metric.impl.ExplicitMetricSpace(distances);
    stream = api.Plex4.createVietorisRipsStream(m_space, max_dimension, max_filtration_value, num_divisions);
    persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);
    intervals = persistence.computeIntervals(stream);
    % barcode0 = edu.stanford.math.plex4.homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, 0);
    barcode1 = edu.stanford.math.plex4.homology.barcodes.BarcodeUtility.getEndpoints(intervals, 1, 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the max_filtration_value
    ck_inf1 = 0;
    if sum(sum(barcode1 > 10^4))
        ck_inf1 = 1;
    end
    count = 1;
    while ck_inf1 == 1
        display(['Increase the max_filtration_value to 3000.']);
        m_space = metric.impl.ExplicitMetricSpace(distances);
        stream = api.Plex4.createVietorisRipsStream(m_space, max_dimension, ...
            max_filtration_value+count*500, num_divisions+count*500);
        persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);
        intervals = persistence.computeIntervals(stream);
        % barcode0 = edu.stanford.math.plex4.homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, 0);
        barcode1 = edu.stanford.math.plex4.homology.barcodes.BarcodeUtility.getEndpoints(intervals, 1, 0);
        
        ck_inf1 = 0;
        if sum(sum(barcode1 > 10^4))
            ck_inf1 = 1;
        end
        count = count + 1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % Change the filtration value.
    [Tree, pred] = graphminspantree(sparse(D));
    
    barcode0 = sort(nonzeros(Tree)); % sval(barcode0(2:end,2));
    barcode1(:,1) = sval(barcode1(:,1));
    barcode1(:,2) = sval(barcode1(:,2));
    
    barcode0_D{g} = barcode0;
    barcode1_D{g} = barcode1;
    display(num2str([g]));
end
clear intervals; clear m_space; 
clear persistence; clear stream; 
``` 

`barcode0_D` and `barcode1_D` have 4 cells, and each cell is a matrix of which row represents the zeroth and first holes, respectively, with the birth and death of thresholds.   
The persistence diagram in `barcode1_D` can be plotted as follows: 

