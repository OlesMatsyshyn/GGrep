 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
update_BS = 0;
%% Band structure computation
disp('Computing effective Hamiltonians')
if update_BS == 1
    H_Effective_T
    %H_Effective_B
end

%% Computing tables for the BdG Hamiltonian
load('HT.mat')
%load('HB.mat')
parameters_minimal

% As we remember
% H_BdG = @(kx,ky,q,DA,DB)[g1*H_eff_K1(kx+q,ky)/scale - mu*s0,        [DA,0;0,DB];
%                         [DA,0;0,DB], -g1*conj(H_eff_K2(-kx,-ky))/scale + mu*s0];
nebands = nbands;
nbands  = 2*nbands;
HT_table = zeros(nebands,nebands,nx,ny,q_size);
HB_table = zeros(nebands,nebands,nx,ny);
H0BDG    = zeros(nbands,nbands,  nx,ny,q_size);
VX_table = zeros(nbands,nbands,  nx,ny,q_size);
disp('Preparing the H and V tables')
for iqx = 1 : q_size
    for i = 1 : nx
        for j = 1 : ny
            HT_table(:,:,i,j,iqx) =  g1*H_eff_T(:,:,i,j,iqx,1)/scale - mu*s0;
            HB_table(:,:,i,j) = -g1*conj(H_eff_T(:,:,end-i+1,end-j+1,fix(length(qx)/2)+1,2))/scale + mu*s0;
            % Regardless of how we include q:
            H0BDG(:,:,i,j,iqx) = [HT_table(:,:,i,j,iqx),  [0,0;0,0];...
                                 [0,0;0,0],      HB_table(:,:,i,j)];
        end
    end
end
for iqx = 3 : q_size-2
    for i = 1 : nx
        for j = 1 : ny
            VX_table(:,:,i,j,iqx) = -H0BDG(:,:,i,j,iqx+2)+8*H0BDG(:,:,i,j,iqx+1)...
                                  -8*H0BDG(:,:,i,j,iqx-1)+  H0BDG(:,:,i,j,iqx-2);
        end
    end
end
VX_table = VX_table / (g1*(qx(2)-qx(1)));


DAq = qx*0+0.035; % seed value
DBq = qx*0+0.015; % seed value
Jxq_v1 = qx*0;
Jxq_v2_m = qx*0;
Jxq_v2_s = qx*0;
Jyq = qx*0;


cprintf('hyper','starting the order parameter canlculation...\n')% XX component (just a nice dislay of the statement)

eps = 10^(-5);
for iqx = 1:q_size
    iter = 0;             % number of iterations from this q completed
    discrA = 1;           % just to inititalize while
    discrB = 1;

    energy = zeros (nbands, nx, ny);
    DAnn   = zeros (nbands, nx, ny);
    DBnn   = zeros (nbands, nx, ny);

    while (discrA > 0.5*eps && abs(DAq(iqx)) > eps)||(discrB > 0.5*eps && abs(DBq(iqx)) > eps)% convergence cryterium
        iter = iter + 1;  % update iteration index
        % Do not count the boundary twise! start from 0 and end at nx-1,
        % or from 1 and end at nx
        for i = 1 : nx % to stay in the limit
            for j = 1 : ny
                [Vec,Val] = eigenshuffle( [HT_table(:,:,i,j,iqx),   [DAq(iqx),0;0,DBq(iqx)]  ;...
                                          [DAq(iqx),0;0,DBq(iqx)],   HB_table(:,:,i,j)    ] );
                energy(:,i,j) = Val(:);    % band, at which kx, at which ky;
                for n = 1 : nbands
                    DAnn(n,i,j)  =  Vec(:,n)'*dDAH_BdG*Vec(:,n);
                    DBnn(n,i,j)  =  Vec(:,n)'*dDBH_BdG*Vec(:,n);
                end
          
            end
        end
        
        DAmem = DAq(iqx);                 % memorise the D before the update
        DBmem = DBq(iqx);
        musk = (abs(energy)<wD);
        % Updating the gap
        DAq(iqx) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DAnn.*musk)));
        DBq(iqx) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DBnn.*musk)));
        %            - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy_K2).*DBnn_K2.*musk_K2)));
         
        discrA = abs((DAq(iqx)-DAmem)); % discrepancy indicator
        discrB = abs((DBq(iqx)-DBmem)); % discrepancy indicator  ', discrepancy = ',num2str(fix(discr*10^6)/10^6)]);
        disp(['iter = ', num2str(iter), ', DA = ',num2str(fix(DAq(iqx)*1000)), ', DB = ',num2str(fix(DBq(iqx)*1000))])
    end
    if DAq(iqx) < 0.5*eps
       DAq(iqx) = 0;
    end
    if DBq(iqx) < 0.5*eps
        DBq(iqx) = 0;
    end
    JxD     = zeros (nbands, nx, ny);
    Jx0     = zeros (nbands, nx, ny);
    energy0 = zeros (nbands, nx, ny);
    for i = 1 : nx
        for j = 1 : ny
            [Vec,Val]   = eigenshuffle( [HT_table(:,:,i,j,iqx),   [DAq(iqx),0;0,DBq(iqx)]  ;...
                                        [DAq(iqx),0;0,DBq(iqx)],   HB_table(:,:,i,j)] );
            [Vecn,Valn] = eigenshuffle( [HT_table(:,:,i,j,iqx), [0,0;0,0]  ;...
                                        [0,0;0,0],     HB_table(:,:,i,j)] );
            energy(:,i,j)  = Val(:);
            energy0(:,i,j) = Valn(:);
            for n = 1 : nbands
                JxD(n,i,j) =  Vec(:,n)'*VX_table(:,:,i,j,iqx)*Vec(:,n);
                Jx0(n,i,j) = Vecn(:,n)'*VX_table(:,:,i,j,iqx)*Vecn(:,n);
            end

        end
    end
    musk  = (abs(energy)<wD);
    Jxq_v1(iqx)   = sum(sum(sum(fFD(energy).*JxD.*musk)));
    Jxq_v2_m(iqx) = sum(sum(sum(fFD(energy).*JxD)));
    Jxq_v2_s(iqx) = sum(sum(sum(fFD(energy0).*Jx0)));
    disp([num2str(fix(1000*(iqx)/q_size)/10),'% finished'])

end
Jxq_v1 = Jxq_v1 * dkxdky; % in e*v0 units
Jxq_v2_m = Jxq_v2_m * dkxdky;
Jxq_v2_s = Jxq_v2_s * dkxdky;
%save('SCDEres0.mat','qx','Jxq_v1','Jxq_v2','wD','mu','U',"DAq","DBq","Delta","DeltaE","Temp",'nx','ny')

figure
plot(qx,DAq,qx,DBq)
xlabel('q_x')
legend('\Delta_A')

figure
plot(qx,real(Jxq_v1),qx,real(Jxq_v2_m-Jxq_v2_s))
legend('J_x')
xlabel('q_x')
ylabel('J_q')



