function Hvals = cumflux_cylinder_inabs_outref(D,tvals,r0,a,b,eiglist)
% calculate the time distribution for particles encountering inner boundary
% at a, assuming they start at r0 and there is a reflecting boundary at b
% eigenvalues for the Sturm-Liouville problem are precalculated in eiglist
% eiglist assumed a row vector
% returns cumulative distribution (probability hit before time t)
% for times in tvals
% D is diffusion coefficient

J0a = besselj(0,eiglist*a);
J1a = besselj(1,eiglist*a);
Y1a = bessely(1,eiglist*a);
J1b = besselj(1,eiglist*b);
Y1b = bessely(1,eiglist*b);
coeff = J0a.^2.*eiglist./(J0a.^2-J1b.^2).*...
    (besselj(0,eiglist*r0).*Y1b - bessely(0,eiglist*r0).*J1b).*...
    (Y1a.*J1b - J1a.*Y1b);

expvals = exp(-eiglist'.^2*D*tvals);
summand = bsxfun(@times,expvals,coeff');
Hvals = 1-pi^2/2*a*sum(summand,1);
end