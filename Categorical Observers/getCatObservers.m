%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function name: getCatObservers
%
% author: Yongmin Park (yp8103@rit.edu)
%
% description: This script generates categorical observers following the
% method proposed by Yuta Asano & Mark Fairchild. Particularly, this script
% utilizes 2019 UN world population data to generate individual observers.
% For further information about categorical observers & individual
% observers, please refer to
% https://www.rit.edu/cos/colorscience/re_AsanoObserverFunctions.php and 
% Asano, Y, Fairchild, MD. Categorical observers for metamerism. Color Res Appl. 2020; 
% 45: 576? 585. https://doi.org/10.1002/col.22493
% 
% parameters: 
%                       n_population: Number of individual observers
%                       fs: FOV (2-deg or 10-deg)
%                       n_cat: Number of categorical observers (must be
%                       smaller than n_population)
% output: 
%                       LMS_cat: Cone fundamentals of the categorical
%                       observers
%                       xyz_CMFs: xyz-like CMFs of the categorical
%                       observers
%                       ages: Ages of the categorical observers
%                       vCat: Physiological parameters of the categorical
%                       observers
% Usage:           [LMS_cat, xyz_CMFs, ages, vCat] = getCatObservers(10000, 2, 10)

function [LMS_cat, xyz_CMFs, ages, vCat] = getCatObservers(n_population, fs, n_cat)
    arguments          
        n_population {mustBeReal}
        fs {mustBeReal}
        n_cat {mustBeReal}
    end

    [LMS_All, var_age, vAll ] = fnc_genMonteCarloObs_UNCensusAgeDist( n_population, fs); % Monte Carlo Observers using US Census Age Distribution
    
    N = length(var_age);             % Number of total observers
    N_cat = n_cat;                               % Number of categorical observers
    LMS_cat = zeros(79, 3, N_cat);
    ages = zeros(N_cat, 1);
    vCat = zeros(N_cat, 8);
    
    % To check if an individual observer is selected as center.
    individual = struct([]);
    
    for i=1:N
        individual(i).ID = i;
        individual(i).center = 0;     % center = 1 (yes) or 0 (no)
    end
    
    % 1st categorical observer: 38 years-old 2-deg observer
    LMS_cat(:, 1, 1) = LMS_All(:, 1, 1);
    LMS_cat(:, 2, 1) = LMS_All(:, 2, 1);
    LMS_cat(:, 3, 1) = LMS_All(:, 3, 1);
    individual(1).center = 1;
    ages(1) = var_age(1);
    vCat(1, :) = vAll(1, :);
    
    % To update clusters at every time.
    cluster = struct([]);
    tcluster = struct([]);            % backup cluster
    
    for i = 1:N_cat
        cluster(i).IDX = 0;            % The index of the center of a given cluster
        cluster(i).objects = [];      % The observers belong to the given cluster
        
        tcluster(i).IDX = 0;
        tcluster(i).objects = [];
    end
    
    cluster(1).IDX = 1;
    tdist = 0;                             % The sum of intra-cluster distances
    mdist = zeros(N_cat, 1);
    
    % To find categorical observers (from 2nd).
    for i = 2: N_cat
        idx = randi([1 N]);     % randomly selcect an observer.
        
        % To avoid selecting an observer who already assigned previously.
        while(individual(idx).center == 1)
            idx = randi([1 N]);
        end
        
        % To fetch the LMS for the selected observer and the observer becomes
        % the center of a given cluster
        LMS_cat(:, :, i) = LMS_All(:, :, idx);
        individual(idx).center = 1;
        cluster(i).IDX = idx;
        
        % step 1. perform an initial clustering.
        for k = 1:N
            % Check if an observer is center or not
            if(individual(k).center == 0)
                p = LMS_All(:, :, k);
                
                for c=1:i
                    q = LMS_cat(:, :, c);
                    mdist(c) = cost_kmedoids(p, q);    % compute the distance between the observer and all the centers in the LMS space.
                end
                
                % This observer belongs to the cluster whose center is the
                % closest from the observer.
                [cost, ic] = min(mdist(1:i));
                
                % Let's update the total intra-cluster distance
                tdist = tdist + cost;
                cluster(ic).objects = [cluster(ic).objects k];
                
            end % end of check if
        end % end of for k
        
        % step 2. find a better center from a given cluster
        for iter = 1:30
            % randomly pick up an observer belonging to the cluster.
            N_cluster = length(cluster(i).objects);
            idx = randi([1 N_cluster]);
            idx = cluster(i).objects(idx);
            
            % Let's do backup the previous center and reset all the observers
            % belonging to the cluster
            for c=1:i
                tcluster(c).IDX = cluster(c).IDX;
                tcluster(c).objects = cluster(c).objects;
                cluster(c).objects = [];
            end
            
            individual(cluster(i).IDX).center = 0;
            individual(idx).center = 1;
            cluster(i).IDX = idx;
            LMS_cat(:, :, i) = LMS_All(:, :, idx);
            
            % The temporal sum of intra-cluster distances
            zdist = 0;
            
            % It's basically the same as the clustering procedure described above (step 1).
            for k = 1:N
                if(individual(k).center == 0)
                    p = LMS_All(:, :, k);
                    
                    for c=1:i
                        q = LMS_cat(:, :, c);
                        mdist(c) = cost_kmedoids(p, q);
                    end
                    
                    [cost, ic] = min(mdist(1:i));
                    
                    zdist = zdist + cost;
                    cluster(ic).objects = [cluster(ic).objects k];
                end
            end
            
            % Let's compute the difference between the previous sum of
            % intra-cluster distances and the present sum of intra-cluster
            % distances
            S = zdist - tdist;
            
            disp([num2str(i) ', ' num2str(iter) ', ' num2str(S)]);
            
            if(S < 0) % If a new clustering gives a smaller sum of intra-cluster distances, then update the sum of intra-cluster distances
                tdist = zdist;
            else % Otherwise, do rollback.
                individual(idx).center = 0;
                for c=1:i
                    cluster(c).IDX = tcluster(c).IDX;
                    cluster(c).objects =  tcluster(c).objects;
                    LMS_cat(:, :, i) = LMS_All(:, :, cluster(c).IDX);
                    individual(cluster(c).IDX).center = 1;
                end
            end % end of else
        end  % end of for iter
    end
    
    % LMS to xyzCMF for 2-deg observer (by Yuta Asano)
    if(fs == 2)        
        M = [0.4151 -0.2424 0.0425;
            0.1355 0.0833 -0.0043;
            -0.0093 0.0125 0.2136];
    elseif(fs == 10)
        M = [0.4499 -0.2630 0.0460;
            0.1617 0.0726 -0.0011;
            -0.0036 0.0054 0.2291];
    else
        M = zeros(3, 3);
        disp('The function only supports either 2 FOV or 10 FOV');
    end
    
    xyz_CMFs = nan(79, 3, N_cat);
    
    for i = 1:N_cat
        t_LMS = LMS_cat(:, :, i);
        xyz_CMFs(:, :, i) = (M * t_LMS')';
    end
    
    % Age distribution
    for i = 1:N_cat
        ages(i) = var_age(cluster(i).IDX);
        vCat(i, :) = vAll(cluster(i).IDX, :);
    end   

end