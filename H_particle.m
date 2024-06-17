function H_eff_T = H_particle(kxi,kyj,qxk,xi)
    % Testings are commented
    parameters_minimal
    % xi = -1;
    % eig_disr = 0;
    band       = zeros (nbands, nbands);
    band_eff   = zeros (nbands, nbands);
    energy     = zeros (nbands, 1);
    energy_eff = zeros (nbands, 1);
    % energy_eff_test = zeros (nbands, 1);
    % Here H are in the units og g1
    H_eff = @(kx,ky)[Delta*(1-kx.^2-ky.^2)      , (xi*kx-1j*ky).^N;  
                       (xi*kx+1j*ky).^N, -Delta*(1-kx.^2-ky.^2)      ]...             % H_0 + H_D
                      + (0.5*g2/g1 - 2 * g3 * (kx.^2+ky.^2) / g0 )*sx...             % H_t
                      + ((delta+D2)/g1 - (2 * g4/g0 + 3*D2/g1) * (kx.^2+ky.^2) )*s0; % H_s + H_s'

    % Define the coupling with the EF, N = 3.
    V_0 = @(kx,ky)lambda^2*[-2*Delta      , 6*(xi*kx-1j*ky);  
                              6*(xi*kx+1j*ky), 2*Delta      ]...% H_0 + H_D
                      - 4*lambda^2 * (g3 / g0 )*sx...            % H_t
                      - 2*lambda^2 * (2 * g4/g0 + 3*D2/g1) * s0; % H_s + H_s'
    
    Vp1 = @(kx,ky) 3*1i*lambda*[0    , (xi*ex-1j*ey)*(xi*kx-1j*ky).^2;  
                                  (xi*ex+1j*ey)*(xi*kx+1j*ky).^2, 0     ]...     % H_0 
                      - 2*1j*lambda*Delta* ( kx*ex+ky*ey )*sz            ...     % H_D
                      - 4*1j*lambda* (g3 / g0 )* ( kx*ex+ky*ey )*sx...           % H_t
                      - 2*1j*lambda* (2 * g4/g0 + 3*D2/g1)* ( kx*ex+ky*ey ) * s0;% H_s + H_s'
    
    Vm1 = @(kx,ky)-3*1i*lambda*[0    , (xi*exc-1j*eyc)*(xi*kx-1j*ky).^2;  
                                  (xi*exc+1j*eyc)*(xi*kx+1j*ky).^2, 0     ]...     % H_0 
                      + 2*1j*lambda*Delta* ( kx*exc+ky*eyc )*sz            ...     % H_D
                      + 4*1j*lambda* (g3 / g0 )* ( kx*exc+ky*eyc )*sx...           % H_t
                      + 2*1j*lambda* (2 * g4/g0 + 3*D2/g1)* ( kx*exc+ky*eyc ) * s0;% H_s + H_s'



    H_eff_T   = zeros (nbands, nbands);



    [Vec,Val] = eigenshuffle(H_eff(kxi+qxk,kyj));
    band(:,:) = Vec;       % column eigenvectors at which kx, at which ky, ;
    energy(:) = Val(:);    % band, at which kx, at which ky;
    % Testings are commented
    % eig_disr = eig_disr + TestEing(H_eff(kxi+qxk,kyj),energy(:),band(:,:));
    % if eig_disr > 10^(-5)
    %     disp('eigenality error')
    % end
    for n = 1 : nbands % a = n, d = m, b = p
       energy_eff(n) = energy(n) + Vec(:,n)'*V_0(kxi,kyj)*Vec(:,n);
       band_eff(:,n) = band(:,n);
       for m = 1 : nbands
            if n == m
                continue;
            end
            c0 = Vec(:,m)'*V_0(kxi,kyj)*Vec(:,n)/(energy(n)-energy(m));
            c1 = - Vec(:,m)'*Vp1(kxi,kyj)*Vec(:,n)*Vec(:,n)'*Vm1(kxi,kyj)*Vec(:,n)...
                 + Vec(:,m)'*Vm1(kxi,kyj)*Vec(:,n)*Vec(:,n)'*Vp1(kxi,kyj)*Vec(:,n);  
            c1 = (1/Omega)*c1/(energy(n)-energy(m));
            band_eff(:,n) = band_eff(:,n) + (c0+c1) * band_eff(:,m);
       end
       band_eff(:,n) = band_eff(:,n)/norm(band_eff(:,n));
    end

    % if sum(sum(sum(sum(imag(energy_eff))))) > 0.001
    %     disp('energy is complex error')
    % end
    % if sum(sum(sum(sum(sum(isnan(band_eff)))))) > 0
    %     disp('eigenstate is nan')
    % end
    
    

    for n = 1 : nbands
        H_eff_T = H_eff_T + energy_eff(n)*band_eff(:,n)*band_eff(:,n)';
    end
    % [Vec,Val] = eigenshuffle(H_eff_T(:,:));
    % energy_eff_test(:) = Val(:);    % band, at which kx, at which ky;
    % eig_disr = eig_disr + TestEing(H_eff_T(:,:),energy_eff_test(:),Vec(:,:));
    % if eig_disr > 10^(-5)
    %    disp('eigenality error')
    % end
end