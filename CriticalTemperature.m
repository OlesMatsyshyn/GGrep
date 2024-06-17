 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
%% Band structure computation
SimpleBZ

sizeT = 20;
Tarr = linspace(0.0000000001,0.05,sizeT);

%DAq = Tarr*0+0.05*g1/scale; % seed value
%DBq = Tarr*0+0.01*g1/scale; % seed value
DAq = Tarr*0+0.035; % seed value
DBq = Tarr*0+0.015; % seed value
Jq  = Tarr*0;

kx = linspace(-1,1,nx);
ky = linspace(-1,1,ny);
[KX,KY] = meshgrid(ky,kx);

dkxdky = ( kx(2)-kx(1) ) * ( ky(2)-ky(1) )/(2*pi)^2;
eps = 10^(-6);
tic
cprintf('hyper','Starting the canlculation... Hold on... It"s not so simple...\n')% XX component (just a nice dislay of the statement)
for it = 1:sizeT
    Temp = Tarr(it);% The system temperature
    disp(['Temperature array index ', num2str(it), ' started'])
    fFD  = @(x) 1/2-tanh(0.5*x/Temp)/2;
    iter = 0;             % number of iterations from this q completed

    discrA = 1;           % just to inititalize while

    energy = zeros (nbands, nx, ny);
    DAnn   = zeros (nbands, nx, ny);
    DBnn   = zeros (nbands, nx, ny);
    % Defining the part of space involved in pairing
    if it > 1
        DAq(it) = DAq(it-1);
        DBq(it) = DBq(it-1);
    end
    while (discrA > 5*eps && abs(DAq(it)) > eps)||(discrB > 5*eps && abs(DBq(it)) > eps)% convergence cryterium
        iter = iter + 1;  % update iteration index
        for i = 1 : nx
            for j = 1 : ny
                [Vec,Val] = eigenshuffle(H_BdG(kx(i),ky(j),0,DAq(it),DBq(it)));
                energy(:,i,j) = Val(:);    % band, at which kx, at which ky;
                for n = 1 : nbands
                    DAnn(n,i,j)  =  Vec(:,n)'*dDAH_BdG*Vec(:,n);
                    DBnn(n,i,j)  =  Vec(:,n)'*dDBH_BdG*Vec(:,n);
                end
                
            end
        end  
                    
        musk = double(abs(energy)<wD);

        DAmem = DAq(it);                 % memorise the D before the update
        DBmem = DBq(it);
        % Updating the gap
        DAq(it) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DAnn.*musk)));
        DBq(it) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DBnn.*musk)));
         
        discrA = abs((DAq(it)-DAmem)); % discrepancy indicator
        discrB = abs((DBq(it)-DBmem)); % discrepancy indicator
        discr = discrA + discrB; 
        disp(['it = ', num2str(iter), ' fin, DA = ', num2str(fix(DAq(it)*1000)),...
              ' DB = ', num2str(fix(DBq(it)*1000)), ' discrepancy = ',num2str((discr*10^3))]);      
     
    end
    if DAq(it) < eps %&& DBq(it) < 10^(-5)
        DAq(it:end) = 0;
        DBq(it:end) = 0;
        break;
    end
end

plot(Tarr,DAq,Tarr,DBq)
xlabel('T')
legend('D_A')
% Comments:
% Convergence with the discr 10^(-2) and 10^(-3) was tested
% The difference is insignificant, but -2 is much faster
