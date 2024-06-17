 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
parameters_minimal

nebands = nbands;
nbands  = 2*nbands;
HT_table = zeros(nebands,nebands,nx,ny,q_size);
Jx_table = zeros(nebands,nebands,nx,ny,q_size);
HB_table = zeros(nebands,nebands,nx,ny);
%VX_table = zeros(nebands,nebands,nx,ny);
disp('Computing the H tables')

for iq = 1 : q_size
    for i = 1 : nx
        for j = 1 : ny
            HT_table(:,:,i,j,iq) =  g1*H_particle(kx(i),ky(j),qx(iq),-1)/scale - mu*s0;
            HB_table(:,:,i,j) =  -g1*conj(H_particle(-kx(i),-ky(j),0,1))/scale + mu*s0;
        end
    end
end

DAq = qx*0+0.035; % seed value
DBq = qx*0+0.015; % seed value
Jxq_v1 = qx*0;
Jxq_v2 = qx*0;

dkxdky = ( kx(2)-kx(1) ) * ( ky(2)-ky(1) )/(2*pi)^2;

cprintf('hyper','starting the order parameter canlculation...\n')% XX component (just a nice dislay of the statement)

eps = 10^(-5);
for k = 1:q_size
    iqx = k;
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
                                          [DAq(iqx),0;0,DBq(iqx)],       HB_table(:,:,i,j)] );
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
    
    disp([num2str(fix(1000*(iqx)/q_size)/10),'% finished'])

end
%save('SCDEres0.mat','qx','Jxq_v1','Jxq_v2','wD','mu','U',"DAq","DBq","Delta","DeltaE","Temp",'nx','ny')
dq = 0.01; 
dxH_particle = @(kx,ky,q) (-H_particle(kx,ky,q+2*dq,-1)+8*H_particle(kx,ky,q+dq,-1)...
                            -8*H_particle(kx,ky,q-  dq,-1)+  H_particle(kx,ky,q-2*dq,-1))/(12*dq*scale);

for iqx = 1 : q_size
    JxD     = zeros (nbands, nx, ny);
    Jx0     = zeros (nbands, nx, ny);
    energy0 = zeros (nbands, nx, ny);
    for i = 1 : nx
        for j = 1 : ny
            [Vec,Val]   = eigenshuffle( [HT_table(:,:,i,j,iqx),   [DAq(iqx),0;0,DBq(iqx)]  ;...
                                        [DAq(iqx),0;0,DBq(iqx)],       HB_table(:,:,i,j)] );
            [Vecn,Valn] = eigenshuffle( [HT_table(:,:,i,j,iqx), [0,0;0,0]  ;...
                                        [0,0;0,0],     HB_table(:,:,i,j)] );
            energy(:,i,j)  = Val(:);
            energy0(:,i,j) = Valn(:);
            for n = 1 : nbands
                JxD(n,i,j) =  Vec(:,n)'*[dxH_particle(kx(i),ky(j),qx(iqx)), zeros(2,2);zeros(2,2),zeros(2,2)]*Vec(:,n);
                Jx0(n,i,j) = Vecn(:,n)'*[dxH_particle(kx(i),ky(j),qx(iqx)), zeros(2,2);zeros(2,2),zeros(2,2)]*Vecn(:,n);
            end

        end
    end
    musk  = (abs(energy)<wD);
    musk0 = (abs(energy0)<wD);
    Jxq_v1(iqx) = sum(sum(sum(fFD(energy).*JxD.*musk)));
    Jxq_v2(iqx) = sum(sum(sum(fFD(energy).*JxD)))-sum(sum(sum(fFD(energy0).*Jx0)));
    disp([num2str(fix(1000*(iqx)/q_size)/10),'% of jx finished'])
end



Jxq_v1 = Jxq_v1 * dkxdky/scale; % in e*v0 units
Jxq_v2 = Jxq_v2 * dkxdky/scale;
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



