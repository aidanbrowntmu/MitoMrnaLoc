function [tvals,Hvals_avg] = getTimeDistrib_refabs_uniform_func(inner_radius,outer_radius)

% clear all;
% close all;
% 
% inner_radius = 0.25;
% outer_radius = 0.35;

disc_number = 100;
fraction_values = linspace(1/(disc_number+1),disc_number/(disc_number+1),disc_number);

for j=1:disc_number
    %j
    fraction = fraction_values(j);
    initial_radius = inner_radius + fraction*(outer_radius - inner_radius);

    [tvals,Hvals] = getTimeDistrib_refabs_func(inner_radius,outer_radius,initial_radius);
    for i=1:length(Hvals)
      if(Hvals(i)<1e-10)
        Hvals(i) = 0;
      end
    end

    radius_keep(j) = initial_radius;
    tvals_keep(j,:) = tvals;
    Hvals_keep(j,:) = Hvals;
end

Hvals_avg(1:length(Hvals)) = 0;
radius_total = 0;
for j=1:disc_number
  Hvals_avg = Hvals_avg + Hvals_keep(j,:)*radius_keep(j);
  radius_total = radius_total + radius_keep(j);
end
Hvals_avg_before = Hvals_avg;

Hvals_avg = Hvals_avg/radius_total;