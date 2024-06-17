function H = H_BdG(kx,ky,qx,DAq,DBq)
    parameters_minimal
    % Here H are in the units og g1
    H_eff = @(xi,kx,ky)[Delta*(1-kx.^2-ky.^2)      , (xi*kx-1j*ky).^N;  
                       (xi*kx+1j*ky).^N, -Delta*(1-kx.^2-ky.^2)      ]...             % H_0 + H_D
                      + (0.5*g2/g1 - 2 * g3 * (kx.^2+ky.^2) / g0 )*sx...             % H_t
                      + ((delta+D2)/g1 - (2 * g4/g0 + 3*D2/g1) * (kx.^2+ky.^2) )*s0; % H_s + H_s'

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



H_eff_T   = zeros (nbands, nbands, q_size, 2, nx, ny);
for iq = 1 : length(qx)
    eig_disr = 0;
    band       = zeros (nbands, nbands, 2, nx, ny);
    band_eff   = zeros (nbands, nbands, 2, nx, ny);
    energy     = zeros (nbands, 2, nx, ny);
    energy_eff = zeros (nbands, 2, nx, ny);
    for xi = -1 : 2 : 1 
        if xi == -1
            id = 1;
        else
            id = 2;
        end
        for i = 1 : nx % to stay in the limit
            for j = 1 : ny
                [Vec,Val] = eigenshuffle(H_eff(xi,kx(i)+qx(iq),ky(j)));
                band(:,:,id,i,j) = Vec;       % column eigenvectors at which kx, at which ky, ;
                                              % eigenvector vi = band(:,i,kx,ky);
                energy(:,id,i,j) = Val(:);    % band, at which kx, at which ky;
                eig_disr = eig_disr + TestEing(H_eff(xi,kx(i)+qx(iq),ky(j)),energy(:,id,i,j),band(:,:,id,i,j));
                if eig_disr > 10^(-5)
                    disp('eigenality error')
                end
                for n = 1 : nbands % a = n, d = m, b = p
                   energy_eff(n,id,i,j) = energy(n,id,i,j) + Vec(:,n)'*V_0(xi,kx(i),ky(j))*Vec(:,n);
                   band_eff(:,n,id,i,j) = band(:,n,id,i,j);
                   for m = 1 : nbands
                        if n == m
                            continue;
                        end
                        c0 = Vec(:,m)'*V_0(xi,kx(i),ky(j))*Vec(:,n)/(energy(n,id,i,j)-energy(m,id,i,j));
                        c1 = - Vec(:,m)'*Vp1(xi,kx(i),ky(j))*Vec(:,n)*Vec(:,n)'*Vm1(xi,kx(i),ky(j))*Vec(:,n)...
                             + Vec(:,m)'*Vm1(xi,kx(i),ky(j))*Vec(:,n)*Vec(:,n)'*Vp1(xi,kx(i),ky(j))*Vec(:,n);  
                        c1 = (1/Omega)*c1/(energy(n,id,i,j)-energy(m,id,i,j));
                        band_eff(:,n,id,i,j) = band_eff(:,n,id,i,j) + (c0+c1) * band_eff(:,m,id,i,j);
                   end
                   band_eff(:,n,id,i,j) = band_eff(:,n,id,i,j)/norm(band_eff(:,n,id,i,j));
                end
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
    for xi = -1 : 2 : 1 
        if xi == -1
            id = 1;
        else
            id = 2;
        end
        for i = 1 : nx % to stay in the limit
            for j = 1 : ny
                for n = 1 : nbands
                    H_eff_T(:,:,iq,id,i,j) = H_eff_T(:,:,iq,id,i,j) + energy_eff(n,id,i,j)*band_eff(:,n,id,i,j)*band_eff(:,n,id,i,j)';
                end
                [Vec,Val] = eigenshuffle(H_eff_T(:,:,iq,id,i,j));
                energy_eff_test(:,id,i,j) = Val(:);    % band, at which kx, at which ky;
                eig_disr = eig_disr + TestEing(H_eff_T(:,:,iq,id,i,j),energy_eff_test(:,id,i,j),Vec(:,:));
                if eig_disr > 10^(-5)
                    disp('eigenality error')
                end
            end
        end

    end
end

end