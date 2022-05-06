function [pouter,tvals,HvalsM,HvalsP] = getTimeDistrib_absabs_func(inner_radius,outer_radius,initial_radius)

%% calculate cumulative distribution of hitting times
% between two absorbing boundaries

% ---------------------
%parameters: MODIFY THESE
a = inner_radius; % radius of inner absorbing cylinder
b = outer_radius; % radius of outer absorbing cylinder
r0 = initial_radius; % starting radius

D = 1; % diffusivity
nmax = 1000;
%max(500,sqrt(1000*(b-a)^2/(r0-a)^2))*2 % number of eigenvalues to use

% min and max time values, logarithmically spaced times will be calculated
tmin = 1e-6;
tmax = 1e3; 
nt = 1e4; % number of time values

%% get eigenvalues
% get the first few with fine-grained search
%n1 = 20;
%[eiglist1,eigfunc1] = getEigsCyl(a,b,10,0.05);
%dx = eiglist1(end)-eiglist1(end-1);
[eiglist,eigfunc] = getEigsCyl(a,b,nmax,[0,0],0.1);
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
[pouter,HvalsP,HvalsM] = cumflux_cylinder_2abs(D,tvals,r0,a,b,eiglist);

% title('conditional arrival time cdf')
% loglog(tvals,HvalsP,tvals,HvalsM)
% legend('outer bound', 'inner bound')
