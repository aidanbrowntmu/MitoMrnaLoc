function [tvals_first,Hvals_first,tvals_CtoB,Hvals_CtoB,tvals_AtoB,Hvals_AtoB, ...
          pouter_startinner,tvals_BtoAC_startinner,Hvals_BtoA_startinner,Hvals_BtoC_startinner, ...
          pouter_startouter,tvals_BtoAC_startouter,Hvals_BtoA_startouter,Hvals_BtoC_startouter, ...
          tvals_AtoBuniform,Hvals_AtoBuniform] = getDiffusionDistributions(epsilon,mt_radius, ...
          radius1,radius2,cellvolume,mitovolumefrac,mrnaparams)

tubelength = mitovolumefrac*(cellvolume)/(pi*mt_radius*mt_radius); %ask Ximena why we add back the 4.2
%tubelength = mitovolumefrac*(cellvolume+4.2)/(pi*mt_radius*mt_radius);
tuberadius = sqrt(cellvolume/(pi*tubelength));
initial_radius = tuberadius - 0.01;
diffusivity = mrnaparams.diffusivity;
      
[tvals_first,Hvals_first] = getTimeDistrib_absref_func(radius2,tuberadius,initial_radius);
Hvals_first = monotonicLowerLimit(1e-6,Hvals_first);
disp('Initial C -> B dist calculated');

[tvals_CtoB,Hvals_CtoB] = getTimeDistrib_absref_func(radius2-epsilon,tuberadius,radius2+epsilon);
Hvals_CtoB = monotonicLowerLimit(1e-6,Hvals_CtoB);
disp('C -> B dist calculated');

[tvals_AtoB,Hvals_AtoB] = getTimeDistrib_refabs_func(mt_radius,radius1+epsilon,radius1-epsilon);
Hvals_AtoB = monotonicLowerLimit(1e-6,Hvals_AtoB);
disp('A -> B dist calculated');

[pouter_startinner,tvals_BtoAC_startinner,Hvals_BtoA_startinner,Hvals_BtoC_startinner] = getTimeDistrib_absabs_func(radius1-epsilon,radius2+epsilon,radius1+epsilon);
Hvals_BtoA_startinner = monotonicLowerLimit(1e-6,Hvals_BtoA_startinner);
Hvals_BtoC_startinner = monotonicLowerLimit(1e-6,Hvals_BtoC_startinner);
disp('B -> A and C start inner dists calculated');

[pouter_startouter,tvals_BtoAC_startouter,Hvals_BtoA_startouter,Hvals_BtoC_startouter] = getTimeDistrib_absabs_func(radius1-epsilon,radius2+epsilon,radius2-epsilon);
Hvals_BtoA_startouter = monotonicLowerLimit(1e-6,Hvals_BtoA_startouter);
Hvals_BtoC_startouter = monotonicLowerLimit(1e-6,Hvals_BtoC_startouter);
disp('B -> A and C start outer dists calculated');

[tvals_AtoBuniform,Hvals_AtoBuniform] = getTimeDistrib_refabs_uniform_func(mt_radius,radius1+epsilon);
Hvals_AtoBuniform = monotonicLowerLimit(1e-6,Hvals_AtoBuniform);
disp('Uniform A -> B dist calculated');

tvals_first = tvals_first/diffusivity;
tvals_CtoB = tvals_CtoB/diffusivity;
tvals_AtoB = tvals_AtoB/diffusivity;
tvals_BtoAC_startinner = tvals_BtoAC_startinner/diffusivity;
tvals_BtoAC_startouter = tvals_BtoAC_startouter/diffusivity;
tvals_AtoBuniform = tvals_AtoBuniform/diffusivity;



end