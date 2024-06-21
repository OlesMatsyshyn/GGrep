clear all
% System initialization

%% Define the strained hoppings and the parameters of the Hamiltonian
% Basic parameters in eV from arxiv 2109.04345
g0    =   3.1;
g1    =   0.38;
g2    = - 0.015;
g3    =   0.29;
g4    =   0.141;
delta = - 0.0105;
D1    =   0.06/2;   % max 100 meV
D2    = - 0.0023;
scale = 0.001; % using meV units 

lambda = 0.03; % perturbation 
Omega  = 0.001; % 100 smaller?
% D1        mu           D_2/D_1
% 0.1       35.6786      11.38
% 0.03      3.66789      7.3
% 0.05      13.0043      11.90
% 0.04      8.33611      10.25 
% 0.06      17.6726      12.97

Delta  = D1/g1;     
wD     = 1*D1/scale;% Debye frequency  1, 0.75, 0.5 was ok. 
DeltaE = 0*5;
mu     = 17.6726;%0.2*g1/scale;     % chemical potential 5/23 or less Band bottom
U      = 60; % used to be 2*g1/scale % Cooper pairs coupling strength. Yanase has 1.5 of gamma_0 for graphene.
% U = 60; Tc = 0.01 meV

eta    = 1;        % Light chirality
N      = 3;

Temp = 0.0000000001;   % The system temperature

% The grid 

nx = 401; % pick a number that divides by 3 (sum devides by 3)
ny = 401; % best when nx = ny
kx = linspace(-1.1,1.1,nx);
ky = linspace(-1.2,1.2,ny);
dkxdky = ( kx(2)-kx(1) ) * ( ky(2)-ky(1) )/(2*pi)^2;
q_max  = 0.015;
q_size = 41;
qx = linspace(-q_max,q_max,q_size); 

fFD  = @(x) 1/2-tanh(0.5*(x)/Temp)/2;
dfFD = @(x) -(0.25).*sech(0.5*(x)/Temp).^2; % We should later devide by Temp


dDAH_BdG= [zeros(2,2)    , [1,0;0,0];
          [1,0;0,0], zeros(2,2)    ];

dDBH_BdG= [zeros(2,2)    , [0,0;0,1];
          [0,0;0,1], zeros(2,2)    ];
% Here H are in the units og g1


nebands = 2;
nbands  = 4;
HT_table = zeros(nebands,nebands,nx,ny,q_size);
HB_table = zeros(nebands,nebands,nx,ny);
H0BDG    = zeros(nbands,nbands,nx,ny,q_size);

disp('Computing the H tables')

for iq = 1 : q_size
    for i = 1 : nx
        for j = 1 : ny
            HT_table(:,:,i,j,iq) =  g1*     H_particle( kx(i), ky(j),qx(iq),-1) /scale - mu*eye(2);
            HB_table(:,:,i,j)    = -g1*conj(H_particle(-kx(i),-ky(j),0     , 1))/scale + mu*eye(2);

            H0BDG(:,:,i,j,iq) = [HT_table(:,:,i,j,iq), [0,0;0,0];...
                                [0,0;0,0],             [0,0;0,0]];
        end
    end
end

VX_table = zeros(nbands,nbands,nx,ny,q_size);
for iqx = 1 : q_size
    for i = 3 : nx-2
        for j = 1 : ny
            VX_table(:,:,i,j,iqx) = -H0BDG(:,:,i+2,j,iqx)+8*H0BDG(:,:,i+1,j,iqx)...
                                  -8*H0BDG(:,:,i-1,j,iqx)+  H0BDG(:,:,i-2,j,iqx);
        end
    end
end
clear H0BDG
VX_table = VX_table/(12*(kx(2)-kx(1)));
cprintf('hyper','starting the canlculation of the current...\n')% XX component (just a nice dislay of the statement)

%
epsp = 10^(-5);
DAq = qx*0+0.025; % seed value
DBq = qx*0+0.005; % seed value
for iqx = 1:q_size
    iter = 0;             % number of iterations from this q completed
    discrA = 1;           % just to inititalize while
    discrB = 1;

    energy = zeros (nbands, nx, ny);
    DAnn   = zeros (nbands, nx, ny);
    DBnn   = zeros (nbands, nx, ny);

    while (discrA > 0.5*epsp && abs(DAq(iqx)) > epsp)||(discrB > 0.5*epsp && abs(DBq(iqx)) > epsp)% convergence cryterium
        iter = iter + 1;  % update iteration index
        % Do not count the boundary twise! start from 0 and end at nx-1,
        % or from 1 and end at nx
        for i = 1 : nx % to stay in the limit
            for j = 1 : ny
                [Vec,Val] = eigenshuffle( [HT_table(:,:,i,j,iqx),   [DAq(iqx),0;0,DBq(iqx)]  ;...
                                          [DAq(iqx),0;0,DBq(iqx)],  HB_table(:,:,i,j)     ] );
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
    if DAq(iqx) < 0.5*epsp
       DAq(iqx) = 0;
    end
    if DBq(iqx) < 0.5*epsp
        DBq(iqx) = 0;
    end
    disp([num2str(fix(1000*(iqx)/q_size)/10),'% finished'])
end

save('DAB.mat','DAq','DBq','qx')
Jxq_v1 = qx*0;
%% == Electric current ==

for iqx = 1 : q_size
    Jx     = zeros (nbands, nx, ny);
    energy = zeros (nbands, nx, ny);
    for i = 1 : nx
        for j = 1 : ny
            [Vec,Val] = eigenshuffle( [HT_table(:,:,i,j,iqx),     [DAq(iqx),0;0,DBq(iqx)]  ;...
                                        [DAq(iqx),0;0,DBq(iqx)],    HB_table(:,:,i,j)     ] );
            energy(:,i,j) = Val(:);
            for n = 1 : nbands
                Jx(n,i,j) = Vec(:,n)'*VX_table(:,:,i,j,iqx)*Vec(:,n);
            end
        end
    end
    musk = (abs(energy)<wD);
    Jxq_v1(iqx)   = sum(fFD(energy).*Jx.*musk,'all');
    disp([num2str(fix(1000*(iqx)/q_size)/10),'% of jx finished'])
end

Jxq_v1 = Jxq_v1 * dkxdky; % in e*v0 units
%save('SCDEres0.mat','qx','Jxq_v1','Jxq_v2','wD','mu','U',"DAq","DBq","Delta","DeltaE","Temp",'nx','ny')
plot(qx,real(Jxq_v1))

