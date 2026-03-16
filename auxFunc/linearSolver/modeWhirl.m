function kappa = modeWhirl(u,v)
%MODEWHIRL Computes the parameter kappa for the orbit. 
% kappa > 0 : FW mode, kappa < 0 : BW mode
% INPUT:    u: vector of first element of the mode
%           v: vector of second element of the mode
npts = length(u);
ru = abs(u);
rv = abs(v);
eta_u = angle(u);
eta_v = angle(v);
kappa = zeros(1,npts);

for ii = 1:npts
     if ru(ii)* rv(ii) < 1e-16
        kappa(ii) = 0;
    else
        T = [ru(ii)*cos(eta_u(ii)), -ru(ii)*sin(eta_u(ii)); rv(ii)*cos(eta_v(ii)), -rv(ii)*sin(eta_v(ii))];

        H = T * T';
        L = eigs(H);
        lambda = sqrt(L);
        kappa(ii) = min(lambda) / max(lambda);

        diff = mod(eta_v(ii) - eta_u(ii),2*pi);
        if diff > 0 && diff < pi
            kappa(ii) = - kappa(ii);
        end
     end
end

end