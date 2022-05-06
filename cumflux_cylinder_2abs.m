function [pouter,HvalsP,HvalsM] = cumflux_cylinder_inabs_outref(D,tvals,r0,a,b,eiglist)
% calculate the time distribution for particles encountering each boundary
% in a cylindrical region with 2 absorbing boundaries
% pouter = probability of encountering outer before inner
% HvalsP = conditional cumulative distribution for outer boundary
%       = probability arrive by time t, given that outer boundary is
%       encountered first
% HvalsM = conditional cumulative distribution for inner boundary 
%       = probability arrive by time t, given that inner boundary is
%       encountered first
%
% inner boundary at a, outer boundary at b
% particles start at r0
% eigenvalues for the Sturm-Liouville problem are precalculated in eiglist
% eiglist assumed a row vector
% returns cumulative distribution (probability hit before time t)
% for times in tvals
% D is diffusion coefficient

% get the splitting probability (prob hit inner first)
pouter = log(r0/a)/log(b/a);


J0a = besselj(0,eiglist*a);
%Y0a = bessely(0,eiglist*a);
J0b = besselj(0,eiglist*b);
Y0b = bessely(0,eiglist*b);
J00 = besselj(0,eiglist*r0);
Y00 = bessely(0,eiglist*r0);
J1a = besselj(1,eiglist*a);
J1b = besselj(1,eiglist*b);
Y1a = bessely(1,eiglist*a);
Y1b = bessely(1,eiglist*b);

% get outer encounter conditional distribution
coeff = J0a.^2.*eiglist./(J0a.^2-J0b.^2).*...
    (J00.*Y0b - Y00.*J0b).*(J1b.*Y0b - Y1b.*J0b);

expvals = exp(-eiglist'.^2*D*tvals);
summand = bsxfun(@times,expvals,coeff');
HvalsP = 1-pi^2*b/(2*pouter)*sum(summand,1);

% get inner encounter conditional distribution
coeff = J0a.^2.*eiglist./(J0a.^2-J0b.^2).*...
    (J00.*Y0b - Y00.*J0b).*(Y1a.*J0b - J1a.*Y0b);
pinner = 1-pouter;

expvals = exp(-eiglist'.^2*D*tvals);
summand = bsxfun(@times,expvals,coeff');
HvalsM = 1-pi^2*a/(2*pinner)*sum(summand,1);

end