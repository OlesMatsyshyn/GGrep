% System initialization

%% Define the strained hoppings and the parameters of the Hamiltonian
% Basic parameters in eV from arxiv 2109.04345
g0    =   3.1;
g1    =   0.38;
g2    = - 0.015;
g3    =   0.29;
g4    =   0.141;
delta = - 0.0105;
D1    =   0.06;   % 100 meV, there is factor 1/2 later
D2    = - 0.0023;
scale = 0.001; % using meV units 

% D1        mu           D_2/D_1
% 0.1       35.6786      11.38
% 0.03      3.66789      7.3
% 0.05      13.0043      11.90
% 0.04      8.33611      10.25 
% 0.06      17.6726      12.97

Delta  = D1/g1;     
wD     = 0.5*D1/scale;% Debye frequency  0.75, 0.5 was ok. 
DeltaE = 0*5;
mu     = 17.6726;%0.2*g1/scale;     % chemical potential 5/23 or less Band bottom
U      = 60; % used to be 2*g1/scale % Cooper pairs coupling strength. Yanase has 1.5 of gamma_0 for graphene.
% U = 60; Tc = 0.01 meV

eta    = 1;        % Light chirality
N      = 3;

Temp = 0.0000000001;   % The system temperature

fFD  = @(x) 1/2-tanh(0.5*(x)/Temp)/2;
dfFD = @(x) -(0.25).*sech(0.5*(x)/Temp).^2; % We should later devide by Temp



s0 = [1,  0 ;
      0 , 1];
sx = [0,  1 ;
      1,  0];
sy = [0 ,-1j;
      1j, 0];
sz = [1,  0;
      0 ,-1];
sp = [0,  1;
      0 , 0];
sm = [0,  0;
      1 , 0];

% Here H are in the units og g1
H_eff_K1 = @(kx,ky)[Delta*(1-kx.^2-ky.^2)/2      , (kx-1j*ky).^N;  
                   (kx+1j*ky).^N, -Delta*(1-kx.^2-ky.^2)/2      ]...
                   + (0.5*g2/g1 - 2 * g3 * (kx.^2+ky.^2) / g0 )*sx...
                   + ((delta+D2)/g1 - (2 * g4/g0 + 3*D2/g1) * (kx.^2+ky.^2) )*s0;
H_eff_K2 = @(kx,ky)[Delta*(1-kx.^2-ky.^2)/2      , (-kx-1j*ky).^N;  
                   (-kx+1j*ky).^N , -Delta*(1-kx.^2-ky.^2)/2    ]  + ...
                   + (0.5*g2/g1 - 2 * g3 * (kx.^2+ky.^2) / g0 )*sx...
                   + ((delta+D2)/g1 - (2 * g4/g0 + 3*D2/g1) * (kx.^2+ky.^2) )*s0;

H_BdG = @(kx,ky,q,DA,DB)[g1*H_eff_K1(kx+q/2,ky)/scale - mu*s0,       [DA,0;0,DB]   ;
                        [DA,0;0,DB], -g1*conj(H_eff_K2(-kx+q/2,-ky))/scale + mu*s0];

% Beware of the minus in the lower part. 
% The reason for it is that dxH(x) = f(x) is a full function already.
% f(-k + q/2) = d(H(-k+q/2)) / d(-k+q/2), so we don't take out the minus
% from the denominator. Regardless of -d(q/2) H(-k+q/2)= dk H(-k+q/2)



dxH_eff_K1 = @(kx,ky)[-Delta*kx      , N*(kx-1j*ky).^(N-1);  
                     N*(kx+1j*ky).^(N-1), Delta*kx      ]...
                   + (- 4 * g3 * kx / g0 )*sx...
                   + (- 2* (2 * g4/g0 + 3*D2/g1)*kx )*s0...
                   + 2*eta*(N-1) * DeltaE *sz *kx.*(kx.^2+ky.^2).^(N-2)...
                   + N*(4*g3/(N*g0))*DeltaE*sz *((kx+1j*ky).^(N-1)+(kx-1j*ky).^(N-1))...
                   - N*(4*Delta/N)*DeltaE*(sm*(kx+1j*ky).^(N-1)+sp*(kx-1j*ky).^(N-1));


dxH_eff_K2 = @(kx,ky)[-Delta*kx      , -N*(-kx-1j*ky).^(N-1);  
                     -N*(-kx+1j*ky).^(N-1) , Delta*kx    ]  + ...
                   + (- 4 * g3 * kx / g0 )*sx...
                   + (- 2* (2 * g4/g0 + 3*D2/g1)*kx )*s0...
                   -  2*eta*(N-1) * DeltaE *sz *kx.*(kx.^2+ky.^2).^(N-2)...
                   + N*(4*g3/(N*g0))*DeltaE*sz *((kx+1j*ky).^(N-1)+(kx-1j*ky).^(N-1))...
                   - N*(4*Delta/N)*DeltaE*(sm*(kx+1j*ky).^(N-1)+sp*(kx-1j*ky).^(N-1));


dxH_BdG = @(kx,ky,q)[g1*dxH_eff_K1(kx+q/2,ky)/scale, 0*s0                        ; 
                     0*s0              , -g1*conj(dxH_eff_K2(-kx+q/2,-ky))/scale]; 

dDAH_BdG= [zeros(2,2)    , [1,0;0,0];
          [1,0;0,0], zeros(2,2)    ];

dDBH_BdG= [zeros(2,2)    , [0,0;0,1];
          [0,0;0,1], zeros(2,2)    ];

nbands = 2*length(H_eff_K1(1,1));