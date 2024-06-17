 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
%% Band structure computation
disp('Initialising')
SimpleBZ


% Ok value 1.1869e-05
% Grid density
nx = 1017; % pick a number that divides by 3 (sum devides by 3)
ny = 1007; % best when nx = ny
kx = linspace(-1.,1.,nx);
ky = linspace(-1.,1.,ny);
l = 1;
kxd = linspace(-1.,1.,l*nx);
kyd = linspace(-1.,1.,l*ny);

%% Computing tables for the BdG Hamiltonian
q_max  = 0.015;
iq_max = fix(0.5*nx*q_max/max(kx)); 
% As we remember
% H_BdG = @(kx,ky,q,DA,DB)[g1*H_eff_K1(kx+q,ky)/scale - mu*s0,        [DA,0;0,DB];
%                         [DA,0;0,DB], -g1*conj(H_eff_K2(-kx,-ky))/scale + mu*s0];
nebands = length(H_eff_K1(1,1));

HT_table = zeros(nebands,nebands,nx,ny);
HB_table = zeros(nebands,nebands,nx,ny);
VX_table = zeros(nebands,nebands,nx,ny);
disp('Computing the H tables')
for i = 1 : nx
    for j = 1 : ny
        HT_table(:,:,i,j) = g1*H_eff_K1(kx(i),ky(j))/scale - mu*s0;
        HB_table(:,:,i,j) = -g1*conj(H_eff_K2(-kx(i),-ky(j)))/scale + mu*s0;
        VX_table(:,:,i,j) = g1*dxH_eff_K1(kx(i),ky(j))/scale;
    end
end
% then effectively Hamiltonian is
%H_BdG = @(i,j,k,DA,DB)[HT_table(i+k,j),   [DA,0;0,DB];
%                      [DA,0;0,DB],     HB_table(i,j)];
%%

qx = (-iq_max:iq_max)*(kx(2)-kx(1));
nqx = length(qx);
DAq = qx*0+0.035; % seed value
DBq = qx*0+0.015; % seed value
Jxq_v1 = qx*0;
Jxq_v2 = qx*0;
Jyq = qx*0;

dkxdky = ( kx(2)-kx(1) ) * ( ky(2)-ky(1) )/(2*pi)^2;

filename = 'progress.txt';
fid = fopen( filename, 'wt' );
fprintf( fid,'%d\n',1);
fclose(fid);

filenameR = 'results.txt';
fidR = fopen( filenameR, 'wt' );
fprintf(fidR,'q DA DB Jx Jy \n\n');
fclose(fidR);


cprintf('hyper','Starting the order parameter canlculation...\n')% XX component (just a nice dislay of the statement)

eps = 10^(-4);
for k = -iq_max:iq_max
    iqx = k+iq_max+1;
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
        for i = 1+iq_max : nx-iq_max % to stay in the limit
            for j = 1 : ny
                [Vec,Val] = eigenshuffle( [HT_table(:,:,i+k,j),   [DAq(iqx),0;0,DBq(iqx)]  ;...
                                          [DAq(iqx),0;0,DBq(iqx)],     HB_table(:,:,i,j)] );
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
    for i = 1+iq_max : nx-iq_max
        for j = 1 : ny
            [Vec,Val]   = eigenshuffle( [HT_table(:,:,i+k,j),   [DAq(iqx),0;0,DBq(iqx)]  ;...
                                        [DAq(iqx),0;0,DBq(iqx)],     HB_table(:,:,i,j)] );
            [Vecn,Valn] = eigenshuffle( [HT_table(:,:,i+k,j),   [0,0;0,0]  ;...
                                        [0,0;0,0],     HB_table(:,:,i,j)] );
            for n = 1 : nbands
                musk = (abs(Val(n))<wD);
                Jxq_v1(iqx) = Jxq_v1(iqx) + fFD(Val(n)).*musk.*Vec(:,n)'*[VX_table(:,:,i+k,j), zeros(2,2);zeros(2,2),zeros(2,2)]*Vec(:,n);
                Jxq_v2(iqx) = Jxq_v2(iqx) + fFD(Val(n)).*Vec(:,n)'*[VX_table(:,:,i+k,j), zeros(2,2);zeros(2,2),zeros(2,2)]*Vec(:,n)...
                                          - fFD(Valn(n)).*Vecn(:,n)'*[VX_table(:,:,i+k,j), zeros(2,2);zeros(2,2),zeros(2,2)]*Vecn(:,n);
            end
        end
    end
     p = update_progress(filename);
     disp([num2str(fix(100*(p)/nqx)),'% finished'])

end
Jxq_v1 = Jxq_v1 * dkxdky/l;
Jxq_v2 = Jxq_v2 * dkxdky/l;
%save('SCDEres0.mat','qx','Jxq_v1','Jxq_v2','wD','mu','U',"DAq","DBq","Delta","DeltaE","Temp",'nx','ny')

figure
plot(qx,DAq,qx,DBq)
xlabel('q_x')
legend('\Delta_A')

figure
plot(qx,real(Jxq_v1),qx,real(Jxq_v2))
legend('J_x')
xlabel('q_x')
ylabel('J_q')



