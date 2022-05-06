function [tvals,Hvals] = getTimeDistrib_absref_func(inner_radius,outer_radius,initial_radius)

%% calculate cumulative distribution of hitting times
% between an outer reflecting and inner absorbing cylinder

% ---------------------
%parameters: MODIFY THESE
a = inner_radius; % radius of inner absorbing cylinder
b = outer_radius; % radius of outer absorbing cylinder
r0 = initial_radius; % starting radius


D = 1; % diffusivity
nmax = 1000;
%round(max(500,sqrt(1000*(b-a)^2/(r0-a)^2))) % number of eigenvalues to use

% min and max time values, logarithmically spaced times will be calculated
tmin = 1e-6;
tmax = 1e3; 
nt = 1e4; % number of time values

%% get eigenvalues
[eiglist,eigfunc] = getEigsCyl(a,b,nmax,[0,1],0.1);
eiglist = eiglist(eiglist>1e-8);

%% check eigenvalues
% xvals = linspace(0,max(eiglist),1000);
% plot(xvals,eigfunc(xvals))
% ylim([-0.01,0.01])
% hold all
% plot(eiglist,zeros(size(eiglist)),'o')
% hold off

%% get cumulative arrival time distribution
tvals = logspace(log10(tmin),log10(tmax),nt);
Hvals = cumflux_cylinder_inabs_outref(D,tvals,r0,a,b,eiglist);

%
%loglog(tvals,Hvals)
