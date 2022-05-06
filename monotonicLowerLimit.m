function distributionnew = monotonicLowerLimit(lowerlimit,distribution)

i = length(distribution);
distributionnew = distribution;
while(i>0)
  if(distribution(i) < lowerlimit)
    distributionnew(1:i) = 0;
  end
  i = i - 1;
end