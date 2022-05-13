function [myMap, com, bin_m, bel_m, exp_m, I, SUM_ENT, ENT_m]  = CRMI( ranges, samscanAngle, Pose, param )
%{
    Copyright (C) 2021  Yang Xu (xuyang94@zju.edu.cn)
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
%}
myResol = param.resol;
% the initial map size in pixels
myMap = zeros(param.size);
mapSize = param.mapSize;
% the origin of the map in pixels
myOrigin = param.origin;
% std dev of z_k
std_dev = param.dev;
maxrange = param.maxrange;
% 4. Log-odd parameters
lo_occ = param.lo_occ;
lo_free = param.lo_free;
lo_max = param.lo_max;
lo_min = param.lo_min;

res_z = param.res_z;

N = size(Pose,2);
numScans = size(samscanAngle);
%% raycast
FLAG = param.flag; % 0 for ijrr map; 1 for wall;
Map = ground_truth_map(myResol, mapSize, FLAG);
%% initialize maps
vec_com = 0.1 * (1:1:10);
com = 0.5 * ones(param.size);
bin_m = zeros(param.size);
bel_m = 1 * ones(10, mapSize * myResol * mapSize * myResol);
exp_m = 0.5 * ones(param.size);
I = zeros(param.size); % CRMI
IND = [];
% map entropy
ENT_m = zeros(param.size);
SUM_ENT = [];

for k = 1:N % for each time with z_k
    disp(['The ', int2str(k),' scan'])
    x = Pose(1, k);
    y = Pose(2, k);
    theta = Pose(3, k);
    
    ix_robot = ceil(x * myResol) + myOrigin(1);
    iy_robot = ceil(y * myResol) + myOrigin(2);
    rays = ranges(:, k);
    %% CRMI initialization
    m_bar = exp_m;
    for i = 1: mapSize * myResol 
    %initializes the MI map with the last map estimate in each loop (entropy in nat)
        I(:,i) = -(exp_m(:,i).*log2(exp_m(:,i)) + (1 - exp_m(:,i)) .* log2(1 - exp_m(:,i))); 
    end
    ENT_m = I;
    %%
    free = [];
    occ = [];
    for i = 1:numScans
        %% find V_ray
        [~, ~, z_ray, z_occ, z_free, ~] = find_nearest_obs(Map, Pose(:,k), samscanAngle(i), param);
        V_ray = z_ray;
        free = [free; z_free];
        occ = [occ; z_occ];
        
        %% CRM
        eta = 0;
        xi = 0.01;
        vec_scm = zeros(length(V_ray));
        for j = 1:length(V_ray)
            com(V_ray(j)) = 2^myMap(V_ray(j))/(1 + 2^myMap(V_ray(j)));
            % update the binary occupancy map
            if com(V_ray(j)) >= 0.7
                bin_m(V_ray(j)) = 1;
            elseif com(V_ray(j)) <= 0.3
                bin_m(V_ray(j)) = -1;
            else
                bin_m(V_ray(j)) = 0;
            end
            % use a simple approximation
            exp_m(V_ray(j)) = (1/10) * sum(vec_com.*bel_m(1:10,V_ray(j))');
        end
        if rays(i) >= maxrange % maxrange case
            prod = 1;
            for j = 1:length(V_ray)
                prod = prod * (1- exp_m(V_ray(j)) / myResol);
                f = prod;
                pz = 1/maxrange;
                % compute SCM p(c_k| z{0:k},x{0:k})
                pc = pz * f;
                vec_scm(j) = pc;
                eta = eta + pc;
            end
            
            vec_scm = insertSort(vec_scm, min(vec_scm));
            eta = eta + 1/maxrange;
        else %  hit case
            for j = 1:length(V_ray)
                if j==1
                    prod = 1;
                else
                    prod = prod * (1- exp_m(V_ray(j-1)) / myResol);
                end
                f = prod * exp_m(V_ray(j));
                [iy,ix] = ind2sub(size(myMap),V_ray(j));
                dist = (1/myResol)*norm([ix_robot,mapSize * myResol + 1 - iy_robot] - [ix,iy]);
                pz = beamSensorModel(dist, rays(i), param);
                pc = pz * f;
                vec_scm(j) = pc;
                eta = eta + pc;
            end
        end
        
        % normalization
        vec_scm = vec_scm / eta;
        %% compute alpha and beta
        for j = 1:length(V_ray)
            beta(j) = 1 + exp_m(V_ray(j))/(1 - exp_m(V_ray(j)))*sum(vec_scm(j + 1:length(V_ray))) - vec_scm(j);
            alpha(j) = (1  - beta(j))/exp_m(V_ray(j));
        end
        
        % compute the map belief
        for j = 1:size(V_ray)
            bel_m(:,V_ray(j)) = bel_m(:,V_ray(j)) .* (alpha(j) * vec_com'  +  beta(j) * ones(10,1));
            % normalization
            bel_m(:,V_ray(j)) = bel_m(:,V_ray(j)) *10 / sum(bel_m(:,V_ray(j)));
        end

        % ray tracing : Find the maxrange cell for each ray
            x_max = maxrange * cos(-samscanAngle(i) + theta) + x ;
            y_max = maxrange * sin(-samscanAngle(i) + theta) + y ;
            ix_max = ceil(x_max * myResol) + myOrigin(1) ;
            iy_max = ceil(y_max * myResol) + myOrigin(2) ;
            ix_max = max(min(ix_max,mapSize * myResol), 1);
            iy_max = max(min(iy_max,mapSize * myResol), 1);
            
            % Find the set of grids where each beam intersects on the map from the source to the maxrange cell: I_ray
            [ix_beam, iy_beam] = mybresenham([ix_robot, iy_robot], [ix_max, iy_max]);
            beamInd = sub2ind(size(myMap), mapSize * myResol + 1 - iy_beam, ix_beam);
            I_ray = beamInd;
                    
        %% update MI map  
            I = update_CRMI(I_ray, rays(i), I, exp_m, m_bar, param);
    end
    
    %% update the log-odds
    
    myMap(free) = myMap(free) - lo_free;
    myMap(occ) = myMap(occ) + lo_occ;
    
    %     Saturate the log-odd values
    myMap = min(myMap,lo_max);
    myMap = max(myMap,lo_min);
    
    imagesc(myMap); hold on;
    plot(ix_robot, mapSize * myResol + 1 - iy_robot, 'rx', 'LineWidth', 3); % indicate robot location
    pause(0.001);
    
    SUM_ENT = [SUM_ENT; sum(sum(ENT_m))];
end

end
function y=insertSort(X, x)
    if x == X(1)
        y = [x; X];
    else
        y = [X; x];
    end
end