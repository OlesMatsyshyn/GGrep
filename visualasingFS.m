load('HT.mat')
load('HB.mat')
parameters_minimal

energy_eff_test = zeros (nbands, 2, nx, ny);
k = 21;
for i = 1 : nx % to stay in the limit
    for j = 1 : ny
        % What are we pairing?
        H1 = H_eff_T(:,:,k,1,i,j);
        H2 = H_eff_T(:,:,fix(q_size/2)+1,2,end-i+1,end-j+1);

        [Vec,Val] = eigenshuffle(H1);
        energy_eff_test(:,1,i,j) = Val(:);    % band, at which kx, at which ky;
        if TestEing(H1,energy_eff_test(:,1,i,j),Vec(:,:)) > 10^(-9)
            assert(false, 'eigenality error T')
        end
        [Vec,Val] = eigenshuffle(H2);
        energy_eff_test(:,2,i,j) = Val(:);    % band, at which kx, at which ky;
        if TestEing(H2,energy_eff_test(:,2,i,j),Vec(:,:)) > 10^(-9)
            assert(false, 'eigenality error B')
        end

    end
end

band_p = zeros (nx, ny);
band_m = zeros (nx, ny);
for n = 1 :nbands
    for i = 1 : nx % to stay in the limit
        for j = 1 : ny
            if abs(g1*energy_eff_test(n,1,i,j)/scale - mu)<wD/10%0.01/2
                band_p(i,j) = 1;
            end
            if abs(g1*energy_eff_test(n,2,i,j)/scale - mu)<wD/10%0.01/2
                band_m(i,j) = 1;
            end
        end
    end
end
figure
imagesc(kx,ky,band_p)
axis equal
figure 
imagesc(kx,ky,band_m)
axis equal

parameters_minimal
figure
plot(kx,squeeze(real(energy_eff_test(1,1,:,fix(end/2)))))
hold on
plot(kx,squeeze(real(energy_eff_test(2,1,:,fix(end/2)))))
plot([min(kx) max(kx)],[mu mu]*scale/g1,'k--')
figure
plot(kx,squeeze(real(energy_eff_test(1,2,:,fix(end/2)))))
hold on
plot(kx,squeeze(real(energy_eff_test(2,2,:,fix(end/2)))))
plot([min(kx) max(kx)],[mu mu]*scale/g1,'k--')
%%
% vp1 = energy_eff_test(1,1,:,:);
% vp2 = energy_eff_test(2,1,:,:);
% vp3 = energy_eff_test(3,1,:,:);
% vp4 = energy_eff_test(4,1,:,:);
% vp5 = energy_eff_test(5,1,:,:);
% vp6 = energy_eff_test(6,1,:,:);
% vp7 = energy_eff_test(5,1,:,:);
% vp8 = energy_eff_test(6,1,:,:);

