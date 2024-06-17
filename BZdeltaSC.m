 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
%% Band structure computation
SimpleBZ
% Good working vertion
nqx = 100;
qx = linspace(-0.015,0.015,nqx);
DAq = qx*0+0.035; % seed value
DBq = qx*0+0.015; % seed value
Jxq_v1 = qx*0;
Jxq_v2 = qx*0;
Jyq = qx*0;
% Ok value 1.1869e-05
% Grid density
nx = 1017; % pick a number that divides by 3 (sum devides by 3)
ny = 1007; % best when nx = ny
kx = linspace(-1.,1.,nx);
ky = linspace(-1.,1.,ny);
l = 2;
kxd = linspace(-1.,1.,l*nx);
kyd = linspace(-1.,1.,l*ny);
dkxdky = ( kx(2)-kx(1) ) * ( ky(2)-ky(1) )/(2*pi)^2;

filename = 'progress.txt';
fid = fopen( filename, 'wt' );
fprintf( fid,'%d\n',1);
fclose(fid);

filenameR = 'results.txt';
fidR = fopen( filenameR, 'wt' );
fprintf(fidR,'q DA DB Jx Jy \n\n');
fclose(fidR);


cprintf('hyper','Starting the canlculation... Hold on... It"s not so simple...\n')% XX component (just a nice dislay of the statement)

eps = 10^(-4);
for iqx = 1:nqx
    iter = 0;             % number of iterations from this q completed
    discrA = 1;           % just to inititalize while
    %discrB = 1;

    energy = zeros (nbands, nx, ny);
    DAnn   = zeros (nbands, nx, ny);
    DBnn   = zeros (nbands, nx, ny);

    while (discrA > 0.5*eps && abs(DAq(iqx)) > eps)||(discrB > 0.5*eps && abs(DBq(iqx)) > eps)% convergence cryterium
        iter = iter + 1;  % update iteration index
        % Do not count the boundary twise! start from 0 and end at nx-1,
        % or from 1 and end at nx
        for i = 1 : nx
            for j = 1 : ny
                [Vec,Val] = eigenshuffle(H_BdG(kx(i),ky(j),qx(iqx),DAq(iqx),DBq(iqx)));
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
    end
    if DAq(iqx) < 0.5*eps
       DAq(iqx) = 0;
    end
    if DBq(iqx) < 0.5*eps
        DBq(iqx) = 0;
    end
    for i = 1 : l*nx
        for j = 1 : l*ny
            [Vec,Val]   = eigenshuffle(H_BdG(kxd(i),kyd(j),qx(iqx),DAq(iqx),DBq(iqx)));
            [Vecn,Valn] = eigenshuffle(H_BdG(kxd(i),kyd(j),qx(iqx),0,0));
            for n = 1 : nbands
                musk = (abs(Val(n))<wD);
                Jxq_v1(iqx) = Jxq_v1(iqx) + fFD(Val(n)).*musk.*Vec(:,n)'*dxH_BdG(kxd(i),kyd(j),qx(iqx))*Vec(:,n);
                Jxq_v2(iqx) = Jxq_v2(iqx) + fFD(Val(n)).*Vec(:,n)'*dxH_BdG(kxd(i),kyd(j),qx(iqx))*Vec(:,n)...
                                          - fFD(Valn(n)).*Vecn(:,n)'*dxH_BdG(kxd(i),kyd(j),qx(iqx))*Vecn(:,n);
                %Jyq(iqx) = Jyq(iqx) + fFD(Val(n)).*Vec(:,n)'*dyH_BdG(k(1),k(2),qx(iqx))*Vec(:,n);
            end
        end
    end
    p = update_progress(filename);
    disp([num2str(fix(100*(p)/nqx)),'% finished'])

    % fileID = fopen(filenameR,'a+');
    % A = [qx(iqx); DAq(iqx); DBq(iqx); Jxq(iqx)* dkxdky/l; Jyq(iqx)* dkxdky/l];
    % fprintf(fileID,'%12.8f %12.8f %12.8 %12.8 %12.8f\n\n',A);
    % fclose(fileID);
end
Jxq_v1 = Jxq_v1 * dkxdky/l;
Jxq_v2 = Jxq_v2 * dkxdky/l;
save('SCDEres0.mat','qx','Jxq_v1','Jxq_v2','wD','mu','U',"DAq","DBq","Delta","DeltaE","Temp",'nx','ny')

figure
plot(qx,DAq,qx,DBq)
xlabel('q_x')
legend('\Delta_A')

figure
plot(qx,real(Jxq_v1),qx,real(Jxq_v2))
legend('J_x')
xlabel('q_x')
ylabel('J_q')



