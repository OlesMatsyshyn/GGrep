
SimpleBZ % initialization

%% Band structure computation
% Electronic Bands        
band_K1   = zeros (nbands/2, nbands/2, nx, ny);
band_K2   = zeros (nbands/2, nbands/2, nx, ny);
energy_K1 = zeros (nbands/2, nx, ny);
energy_K2 = zeros (nbands/2, nx, ny);

kx = linspace(-1,1,nx);
ky = linspace(-1,1,ny);
[KX,KY] = meshgrid(ky,kx);
q = 0.00;
for i = 1 : nx
    for j = 1 : ny
        [Vec,Val] = eigenshuffle((g1/scale)*H_eff_K1(kx(i)+q/2,ky(j)));
        band_K1(:,:,i,j) = Vec;       % column eigenvectors at which kx, at which ky, ;
        energy_K1(:,i,j) = Val(:);    % band, at which kx, at which ky;
        [Vec,Val] = eigenshuffle((g1/scale)*conj(H_eff_K2(-kx(i)+q/2,-ky(j))));
        band_K2(:,:,i,j) = Vec;       % column eigenvectors at which kx, at which ky, ;
        energy_K2(:,i,j) = Val(:);    % band, at which kx, at which ky;
    end
end

unnormalizationfactorBand = zeros(1,nbands);
unorthogonalityfactorBand = zeros(1,nbands);
uneigenalityfctorBand = zeros(1,nbands);

for n = 1 : nbands/2
    for i = 1 : nx
        for j = 1 : ny
            unnormalizationfactorBand(n) = unnormalizationfactorBand(n) + abs(1 - band_K1(:,n,i,j)'*band_K1(:,n,i,j))+ abs(1 - band_K2(:,n,i,j)'*band_K2(:,n,i,j));
            for m = 1 : nbands/2
                if m == n
                    continue;
                end
                unorthogonalityfactorBand(n) = unorthogonalityfactorBand(n) + abs(band_K1(:,n,i,j)'*band_K1(:,m,i,j)) + abs(band_K2(:,n,i,j)'*band_K2(:,m,i,j));
            end    
            uneigenalityfctorBand(n) = uneigenalityfctorBand(n) + sum(abs((g1/scale)*H_eff_K1(kx(i)+q/2,ky(j))*band_K1(:,n,i,j) - energy_K1(n,i,j)*band_K1(:,n,i,j)))...
                +sum(abs(conj((g1/scale)*H_eff_K2(-kx(i)+q/2,-ky(j)))*band_K2(:,n,i,j) - energy_K2(n,i,j)*band_K2(:,n,i,j)));
        end
    end
end
if (sum(unnormalizationfactorBand) + sum(unorthogonalityfactorBand) + sum(uneigenalityfctorBand)> 10^(-5))
    cprintf('err','Too large numerical error...\n');
    return
end

surf(KX,KY,reshape(energy_K1(1,:,:),nx,ny),'EdgeColor','none')
hold on
surf(KX,KY,reshape(energy_K1(2,:,:),nx,ny),'EdgeColor','none')
surf(KX,KY,mu+0*reshape(energy_K1(1,:,:),nx,ny),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0])
surf(KX,KY,mu+wD+0*reshape(energy_K1(1,:,:),nx,ny),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0.5])
surf(KX,KY,mu-wD+0*reshape(energy_K1(1,:,:),nx,ny),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0.5])


%zlim([-0.04 0.04])
% axis equal
xlabel('$a\cdot p_x$','Interpreter','latex')
ylabel('$a\cdot p_y$','Interpreter','latex')
figure

surf(KX,KY,reshape(energy_K2(1,:,:),nx,ny),'EdgeColor','none')
hold on
surf(KX,KY,reshape(energy_K2(2,:,:),nx,ny),'EdgeColor','none')
surf(KX,KY,mu+0*reshape(energy_K2(1,:,:),nx,ny),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0])
surf(KX,KY,mu+wD+0*reshape(energy_K2(1,:,:),nx,ny),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0.5])
surf(KX,KY,mu-wD+0*reshape(energy_K2(1,:,:),nx,ny),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0.5])


%zlim([-0.04 0.04])
% axis equal
xlabel('$a\cdot p_x$','Interpreter','latex')
ylabel('$a\cdot p_y$','Interpreter','latex')

figure
energyslice_K1 = zeros(nx,ny,2);
energyslice_K2 = zeros(nx,ny,2);
for i = 1 : 2
    enl = energy_K1(i,:, :);
    enl(enl<mu-wD) = NaN;
    enl(enl>mu+wD) = NaN;
    energyslice_K1(:,:,i) = double(~isnan(enl));
    enr = energy_K2(i,:, :);
    enr(enr<mu-wD) = NaN;
    enr(enr>mu+wD) = NaN;
    energyslice_K2(:,:,i) = double(~isnan(enr));
end
imagesc(kx,ky,sum(energyslice_K1,3)+sum(energyslice_K2,3))
axis equal
