clear all
parameters_minimal
% Here H are in the units og g1
H_eff = @(xi,kx,ky)[Delta*(1-kx.^2-ky.^2)      , (xi*kx-1j*ky).^N;  
                    (xi*kx+1j*ky).^N, -Delta*(1-kx.^2-ky.^2)      ]...             % H_0 + H_D
                    + (0.5*g2/g1 - 2 * g3 * (kx.^2+ky.^2) / g0 )*sx...             % H_t
                    + ((delta+D2)/g1 - (2 * g4/g0 + 3*D2/g1) * (kx.^2+ky.^2) )*s0; % H_s + H_s'
% The electric field
ex  =  1  / sqrt(2);
ey  =  1i / sqrt(2);
exc =  1  / sqrt(2);
eyc = -1i / sqrt(2);
% Define the coupling with the EF, N = 3.
V_0 = @(xi,kx,ky)lambda^2*[-2*Delta      , 6*(xi*kx-1j*ky);  
                          6*(xi*kx+1j*ky), 2*Delta      ]...% H_0 + H_D
                  - 4*lambda^2 * (g3 / g0 )*sx...            % H_t
                  - 2*lambda^2 * (2 * g4/g0 + 3*D2/g1) * s0; % H_s + H_s'

Vp1 = @(xi,kx,ky) 3*1i*lambda*[0    , (xi*ex-1j*ey)*(xi*kx-1j*ky).^2;  
                              (xi*ex+1j*ey)*(xi*kx+1j*ky).^2, 0     ]...     % H_0 
                  - 2*1j*lambda*Delta* ( kx*ex+ky*ey )*sz            ...     % H_D
                  - 4*1j*lambda* (g3 / g0 )* ( kx*ex+ky*ey )*sx...           % H_t
                  - 2*1j*lambda* (2 * g4/g0 + 3*D2/g1)* ( kx*ex+ky*ey ) * s0;% H_s + H_s'

Vm1 = @(xi,kx,ky)-3*1i*lambda*[0    , (xi*exc-1j*eyc)*(xi*kx-1j*ky).^2;  
                              (xi*exc+1j*eyc)*(xi*kx+1j*ky).^2, 0     ]...     % H_0 
                  + 2*1j*lambda*Delta* ( kx*exc+ky*eyc )*sz            ...     % H_D
                  + 4*1j*lambda* (g3 / g0 )* ( kx*exc+ky*eyc )*sx...           % H_t
                  + 2*1j*lambda* (2 * g4/g0 + 3*D2/g1)* ( kx*exc+ky*eyc ) * s0;% H_s + H_s'



H_eff_T   = zeros (nbands, nbands, nx, ny, q_size,2);
for id = 1 : 2
    if id == 1
        xi = -1;
    elseif id == 2
        xi = 1;
    end
    for iq = 1 : length(qx)
        eig_disr = 0;
        band       = zeros (nbands, nbands, nx, ny);
        band_eff   = zeros (nbands, nbands, nx, ny);
        energy     = zeros (nbands, nx, ny);
        energy_eff = zeros (nbands, nx, ny);
    
        for i = 1 : nx % to stay in the limit
            for j = 1 : ny
                [Vec,Val] = eigenshuffle(H_eff(xi,kx(i)+qx(iq),ky(j)));
                band(:,:,i,j) = Vec;       % column eigenvectors at which kx, at which ky, ;
                                              % eigenvector vi = band(:,i,kx,ky);
                energy(:,i,j) = Val(:);    % band, at which kx, at which ky;
                eig_disr = eig_disr + TestEing(H_eff(xi,kx(i)+qx(iq),ky(j)),energy(:,i,j),band(:,:,i,j));
                if eig_disr > 10^(-5)
                    disp('eigenality error')
                end
                for n = 1 : nbands % a = n, d = m, b = p
                   energy_eff(n,i,j) = energy(n,i,j) + Vec(:,n)'*V_0(xi,kx(i),ky(j))*Vec(:,n);
                   band_eff(:,n,i,j) = band(:,n,i,j);
                   for m = 1 : nbands
                        if n == m
                            continue;
                        end
                        c0 = Vec(:,m)'*V_0(xi,kx(i),ky(j))*Vec(:,n)/(energy(n,i,j)-energy(m,i,j));
                        c1 = - Vec(:,m)'*Vp1(xi,kx(i),ky(j))*Vec(:,n)*Vec(:,n)'*Vm1(xi,kx(i),ky(j))*Vec(:,n)...
                             + Vec(:,m)'*Vm1(xi,kx(i),ky(j))*Vec(:,n)*Vec(:,n)'*Vp1(xi,kx(i),ky(j))*Vec(:,n);  
                        c1 = (1/Omega)*c1/(energy(n,i,j)-energy(m,i,j));
                        band_eff(:,n,i,j) = band_eff(:,n,i,j) + (c0+c1) * band_eff(:,m,i,j);
                   end
                   band_eff(:,n,i,j) = band_eff(:,n,i,j)/norm(band_eff(:,n,i,j));
                end
            end  
    
        end
        if sum(sum(sum(sum(imag(energy_eff))))) > 0.001
            disp('energy is complex error')
        end
        if sum(sum(sum(sum(sum(isnan(band_eff)))))) > 0
            disp('eigenstate is nan')
        end
        
        energy_eff_test = zeros (nbands, 2, nx, ny);
    
    
        for i = 1 : nx % to stay in the limit
            for j = 1 : ny
                for n = 1 : nbands
                    H_eff_T(:,:,i,j,iq,id) = H_eff_T(:,:,i,j,iq,id) + energy_eff(n,i,j)*band_eff(:,n,i,j)*band_eff(:,n,i,j)';
                end
                [Vec,Val] = eigenshuffle(H_eff_T(:,:,i,j,iq,id));
                energy_eff_test(:,i,j) = Val(:);    % band, at which kx, at which ky;
                eig_disr = eig_disr + TestEing(H_eff_T(:,:,i,j,iq,id),energy_eff_test(:,i,j),Vec(:,:));
                if eig_disr > 10^(-5)
                    disp('eigenality error')
                end
            end
        end
    end
end
id = 1;
J_x       = zeros (nbands, nbands, nx, ny, q_size);
for iq = 1 : length(qx)
    for i = 3 : nx-2
        for j = 1 : ny
            J_x(:,:,i,j,iq) = -H_eff_T(:,:,i+2,j,iq,id)+8*H_eff_T(:,:,i+1,j,iq,id)...
                            -8*H_eff_T(:,:,i-1,j,iq,id)+  H_eff_T(:,:,i-2,j,iq,id);
        end
    end
end
J_x = J_x/(12*(kx(2)-kx(1)));
% figure
% hold on
% plot(kx, squeeze(real(energy_eff(1,1,:,fix(end/2)))),'b-')
% plot(kx, squeeze(real(energy_eff(2,1,:,fix(end/2)))),'b-')
% plot(kx, squeeze(real(energy_eff(1,2,end:-1:1,fix(end/2)))),'r-')
% plot(kx, squeeze(real(energy_eff(2,2,end:-1:1,fix(end/2)))),'r-')
save('HT.mat','H_eff_T','J_x')
%clear all