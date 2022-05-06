clear all
close all

%cell geometry parameters
mitovolumefrac = 0.04;
mvf_factor = 0.8;
mitovolumefrac = mitovolumefrac/mvf_factor;
mt_radius = 0.35;
radius1 = mt_radius + 0.025;
radius2 = mt_radius + 0.25;
cellvolume = mvf_factor*42;
epsilon = 0.01;

%general mrna and trajectory parameters
lifetime = 10*60;
mrnaparams.lifetime1 = lifetime;
mrnaparams.lifetime2 = lifetime;
mrnaparams.lifetime3 = lifetime;
nucexitmean = 60;
nucexitsd = 30;
mrnaparams.MTSstart = 100;
mrnaparams.MTSfold = 0.0275;
mrnaparams.diffusivity = 0.1;

%gene-specific parameters
riboparams.elongrate = 5;
riboparams.initrate = 0.1;
mrnaparams.length = 400;

%determines distributions of times to leave zones 1, 2, 3
[tvals_first,Hvals_first,tvals_CtoB,Hvals_CtoB,tvals_AtoB,Hvals_AtoB, ...
                 pouter_startinner,tvals_BtoAC_startinner,Hvals_BtoA_startinner,Hvals_BtoC_startinner, ...
                 pouter_startouter,tvals_BtoAC_startouter,Hvals_BtoA_startouter,Hvals_BtoC_startouter, ...
                 tvals_AtoBuniform,Hvals_AtoBuniform] = ...
                 getDiffusionDistributions(epsilon,mt_radius, ...
                 radius1,radius2,cellvolume,mitovolumefrac,mrnaparams);

%mrna is initially empty of ribosomes
initialocc(1:mrnaparams.length) = 0;

%generate one trajectory, from mRNA synthesis to decay
[time,timerec_z,zonerec,nonstickyfrac,stickytime,stickysoontime]=...
mrnaTimedFold(riboparams,mrnaparams,initialocc,...
tvals_first,Hvals_first,tvals_CtoB,Hvals_CtoB,tvals_AtoB,Hvals_AtoB, ...
      pouter_startinner,tvals_BtoAC_startinner,Hvals_BtoA_startinner,Hvals_BtoC_startinner, ...
      pouter_startouter,tvals_BtoAC_startouter,Hvals_BtoA_startouter,Hvals_BtoC_startouter, ...
      tvals_AtoBuniform,Hvals_AtoBuniform);

%select how long mRNA will remain in nucleus
nucleustime = normrnd(nucexitmean,nucexitsd);
while(nucleustime < 0)
    nucleustime = normrnd(nucexitmean,nucexitsd);
end
  
%this is the fraction of time binding competent, from Fig 1E
stickyfrac = stickytime/(time+nucleustime);

%determine fraction of time spent in each zone during trajectory
timesum_z(1:3) = 0;
for i=1:(length(timerec_z)/2)
    timesum_z(zonerec(2*i)) = timesum_z(zonerec(2*i)) + timerec_z(2*i) - timerec_z(2*i-1);
end
frac(1) = timesum_z(1)/(sum(timesum_z) + nucleustime);
frac(2) = timesum_z(2)/(sum(timesum_z) + nucleustime);
frac(3) = timesum_z(3)/(sum(timesum_z) + nucleustime);

%localized fraction if zones 1 and 2 counted as localized
frac_12 = frac(1) + frac(2);

%localized fraction if binding-competent in zone 1 counted as localized
frac_stickyin1 = time*(timesum_z(1)/sum(timesum_z) - nonstickyfrac)/(time + nucleustime);

%localized fraction if cycloheximide is added
frac_chx = stickysoontime/(time + nucleustime);