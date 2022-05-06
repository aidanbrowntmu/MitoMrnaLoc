function [time,timerec_z,zonerec,nonstickyfrac,stickytime,stickysoontime]=...
    mrnaTimedFold(riboparams,mrnaparams,initialocc,...
  tvals_first,Hvals_first,tvals_CtoB,Hvals_CtoB,tvals_AtoB,Hvals_AtoB, ...
          pouter_startinner,tvals_BtoAC_startinner,Hvals_BtoA_startinner,Hvals_BtoC_startinner, ...
          pouter_startouter,tvals_BtoAC_startouter,Hvals_BtoA_startouter,Hvals_BtoC_startouter, ...
          tvals_AtoBuniform,Hvals_AtoBuniform)
      
riboocc=initialocc; %initial ribosome occupation of mrna
ki_use=riboparams.initrate; %initiation rate
ke_use=riboparams.elongrate; %elongation rate

lastzonetime = 0; %initializing time left last zone
eventcount = 0; %counting events in trajectory
zswitchcount = 0; %counting number of time switched zones
time = 0; %initial time
timelastinterval = 0; %initializing time before current event
degraded = 0; %whether or not mrna has been degraded
zone = 3; %starting zone
needzonepick = 1; %whether need to determine zone leaving time and direction
ribosomesterminated = 0; %initializing count of ribosomes that have finished translation
nonstickytime(1:3) = 0; %initializing nonsticky time spent in each zone
nonstickysoontime(1:3) = 0; %initializing time almost sticky spend in each zone
stickylastinterval = 0; %whether binding competent in last interval

MTSoff = 0; %counting ribosomes that have ended translation
foldingon = 0; %counting ribosomes that have translated the mts region
mts_stick = 0; %number of binding-competent MTSs

timerec_MTS = []; %initializing record of when MTSs become binding competent
MTSstart = mrnaparams.MTSstart; %number of codons for MTS

while (degraded == 0) %end when mrna is degraded
    
    %determining if/when ribosome has binding competent MTS
    nextstickytime = 1e20;
    mts_stick = 0;
    if(foldingon > MTSoff)
        timerec_MTS_present = timerec_MTS(MTSoff+1:foldingon);
        ribosticktime = min(timerec_MTS_present);
        
        if(ribosticktime <= time)
            mts_stick = length(find(timerec_MTS_present<time));
        else
            mts_stick = 0;
            nextstickytime = ribosticktime;
        end
    else
        ribosticktime = 1e20;
    end
    uniform = 0;
    
    %deciding whether its needed to find zone time to leave and direction
    if((zone==1) && (mts_stick>0))
        needzonepick = 0;
        zonetime = 1e20;
    elseif((zone==1) && (zonetime>1e19))
        needzonepick = 1;
        uniform = 1;
    end
    
    %counting time not binding competent
    if((eventcount > 0) && (stickylastinterval < 1))
        nonstickytime(zonelastinterval) = nonstickytime(zonelastinterval) + (time - timelastinterval);
    end
    
    %counting time not past MTS translation part of mrna
    if((eventcount > 0) && (stickysoonlastinterval < 1))
        nonstickysoontime(zonelastinterval) = nonstickysoontime(zonelastinterval) + (time - timelastinterval);
    end
    
    zonelastinterval = zone;
    timelastinterval = time;
    stickylastinterval = mts_stick;
    
    if(foldingon>MTSoff)
      stickysoonlastinterval = 1;
    else
      stickysoonlastinterval = 0;
    end
    
    %determining time to leave zone, direction to leave
    if(needzonepick == 1)
        if(eventcount < 1)
            tvals = tvals_first;
            Hvals = Hvals_first;
        elseif(uniform == 1)
            tvals = tvals_AtoBuniform;
            Hvals = Hvals_AtoBuniform;
        elseif(zone == 1)
            tvals = tvals_AtoB;
            Hvals = Hvals_AtoB;
        elseif(zone == 2)
            if(lastzone == 1)
                if(rand() < pouter_startinner)
                    tvals = tvals_BtoAC_startinner;
                    Hvals = Hvals_BtoC_startinner;
                    newzone = 3;
                else
                    tvals = tvals_BtoAC_startinner;
                    Hvals = Hvals_BtoA_startinner;
                    newzone = 1;
                end
            else
                if(rand() < pouter_startouter)
                    tvals = tvals_BtoAC_startouter;
                    Hvals = Hvals_BtoC_startouter;
                    newzone = 3;
                else
                    tvals = tvals_BtoAC_startouter;
                    Hvals = Hvals_BtoA_startouter;
                    newzone = 1;
                end
            end
        elseif(zone == 3)
            tvals = tvals_CtoB;
            Hvals = Hvals_CtoB;
        end
        
        rannum = rand();
        index = find(Hvals>rannum,1);
        between = (rannum - Hvals(index-1))/(Hvals(index) - Hvals(index-1));
        tau = tvals(index-1) + between*(tvals(index)-tvals(index-1));
        zonetime = time + tau;
        
        needzonepick = 0;
    end
    
    % elongation rates for all ribosomes except last site
    for i=1:(mrnaparams.length-1)
        if((riboocc(i)==1) && (riboocc(i+1)==0))
            ke(i)=ke_use;
        else
            ke(i)=0;
        end
    end
    
    % elongation rate at last site
    if(riboocc(mrnaparams.length)==1)
        ke(mrnaparams.length)=ke_use;
    else
        ke(mrnaparams.length)=0;
    end
    
    %initiation rate
    if riboocc(1)==0
        ki=ki_use;
    else
        ki=0;
    end
    
    %degradation rate
    if(zone == 1)
        kd = 1/mrnaparams.lifetime1;
    elseif(zone == 2)
        kd = 1/mrnaparams.lifetime2;
    else
        kd = 1/mrnaparams.lifetime3;
    end

    %determining gillespie rates
    ktot = sum(ke) + ki + kd;
    delta_t = (1/(ktot))*log(1/rand());
    time = time + delta_t;
    
    %do gillespie if not switching zones or becoming binding competent
    if((time < zonetime) && (time < nextstickytime))
        rannum=rand();
        chosen=0;
        
        beforeoccMTS = riboocc(mrnaparams.MTSstart);
        beforeoccnend = riboocc(length(riboocc));
        
        %choosing to move a ribosome forward one base pair
        if((sum(ke)/ktot) > rannum)
            for kei=1:mrnaparams.length
                if(chosen == 0)
                    if((sum(ke(1:kei))/ktot) > rannum)
                        if(kei < mrnaparams.length)
                            riboocc(kei) = 0;
                            riboocc(kei+1) = 1;
                        else
                            riboocc(kei) = 0;
                            ribosomesterminated = ribosomesterminated + 1;
                        end
                        chosen = 1;
                    end
                end
            end
        end
        afteroccMTS = riboocc(MTSstart);
        afteroccend = riboocc(length(riboocc));
        
        %pick time to become binding competent if just finished MTS
        if(beforeoccMTS == 0) && (afteroccMTS == 1)
            timetostick = 1/mrnaparams.MTSfold * log(1/rand());
            
            foldingon = foldingon + 1;
            timerec_MTS(foldingon) = time + timetostick;
            % list of times that mRNA would be sticky
        end
        
        
        if (beforeoccnend == 1) && (afteroccend == 0)
            MTSoff = MTSoff + 1;
        end
        
        
        %choosing to initiate translation
        if(chosen == 0)
            cumulative = sum(ke) + ki;
            if((cumulative/ktot) > rannum)
                riboocc(1) = 1;
                chosen = 1;
            end
        end
        
        %choosing to degrade mrna
        if(chosen == 0)
            cumulative = cumulative + kd;
            if((cumulative/ktot) > rannum)
                degraded = 1;
                chosen = 1;
                
                timerec_z(2*(zswitchcount+1) - 1) = lastzonetime;
                timerec_z(2*(zswitchcount+1)) = time;
                zonerec(2*(zswitchcount+1)-1) = zone;
                zonerec(2*(zswitchcount+1)) = zone;
                foldrec(2*(zswitchcount+1) - 1) = foldingon - MTSoff;
                foldrec(2*(zswitchcount+1)) = foldingon - MTSoff;
                zswitchcount = zswitchcount + 1;
                
                if(stickylastinterval < 1)
                  nonstickytime(zonelastinterval) = nonstickytime(zonelastinterval) + (time - timelastinterval);
                end
            end
        end
    elseif(zonetime < nextstickytime) %switch zones
        lasttime = lastzonetime;
        lastzonetime = zonetime;
        time = zonetime;
        needzonepick = 1;
        
        if(zone == 1)
            zone = 2;
            lastzone = 1;
        elseif(zone == 2)
            zone = newzone;
            lastzone = 2;
        else
            zone = 2;
            lastzone = 3;
        end
        
        timerec_z(2*(zswitchcount+1) - 1) = lasttime;
        timerec_z(2*(zswitchcount+1)) = time;
        zonerec(2*(zswitchcount+1)-1) = lastzone;
        zonerec(2*(zswitchcount+1)) = lastzone;
        foldrec(2*(zswitchcount+1) - 1) = foldingon - MTSoff;
        foldrec(2*(zswitchcount+1)) = foldingon - MTSoff;
        zswitchcount = zswitchcount + 1;
    else %become binding competent
        time = nextstickytime;
    end
    
    eventcount = eventcount + 1;
    
end

nonstickyfrac = nonstickytime/time;
stickytime = time - sum(nonstickytime);
stickysoontime = time - sum(nonstickysoontime);

end

