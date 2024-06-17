function ret = TestEing(H,energy,psi)
unnormalizationfactorBand = zeros(1,length(H));
unorthogonalityfactorBand = zeros(1,length(H));
uneigenalityfctorBand     = zeros(1,length(H));
for n = 1 : length(H)
    unnormalizationfactorBand(n) = unnormalizationfactorBand(n) + abs(1 - psi(:,n)'*psi(:,n));
     for m = 1 : length(H)
        if m == n
            continue;
        end
        unorthogonalityfactorBand(n) = unorthogonalityfactorBand(n) + abs(psi(:,n)'*psi(:,m));
    end    
    uneigenalityfctorBand(n) = uneigenalityfctorBand(n) + sum(abs(H*psi(:,n) - energy(n)*psi(:,n)));
end
ret = sum(unnormalizationfactorBand+unorthogonalityfactorBand+uneigenalityfctorBand);
end