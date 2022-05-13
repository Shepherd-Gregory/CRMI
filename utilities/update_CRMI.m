function I = update_CRMI(I_ray, ray, I, exp_m, m_bar, param)
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

    for j = 1: size(I_ray)
        h_bar = 0; % Initialize the conditional entropy of map points / grids
        z = 1 / param.res_z; % Initialize distance update variable
        while z<= min(ray, param.maxrange) % Calculate the MI of each cell and current beam
            % Calculate marginalprobability prob z
            % p（z|M=0）Calculated by uniform distribution
            p1 = param.resol/param.maxrange; 
            p2 = 0;
            for kk = 1: size(I_ray)
                p1 = p1 * (1 - exp_m(I_ray(kk)));
                if kk==1
                    prod2 = 1;
                else
                    prod2 = prod2 * (1- exp_m(I_ray(kk-1))/ param.resol);
                end
                f2 = prod2 * exp_m(I_ray(kk));
                % p(z_k|m_i,x_k)
                pz2 = beamSensorModel(z, min(ray, param.maxrange), param);
                p2 = p2 + pz2 * f2;
            end
            xi = p2 + p1;
            if j==1
                prod3 = 1;
            else
                prod3 = prod3 * (1- exp_m(I_ray(j-1))/ param.resol);
            end
            f3 = prod3 * exp_m(I_ray(j));
            % p(z_k|m_i,x_k)
            pz3 = beamSensorModel(z, min(ray, param.maxrange), param);
            m_bar(I_ray(j)) = pz3 * f3 / xi;
            dh_bar = xi * (m_bar(I_ray(j)) * log2(m_bar(I_ray(j))) + (1 - m_bar(I_ray(j))) * log2(1 - m_bar(I_ray(j))));
            h_bar = h_bar + dh_bar;
            z = z + 1 / param.res_z;
        end
        I(I_ray(j)) = I(I_ray(j)) + h_bar / param.res_z ;
        I(I_ray(j)) = max(0.01, I(I_ray(j)));        
    end
end