inp=arn2_timetest.inp
sysconf
rendld=15
airyfl=f
jmax=60
spac=0.1
energ=500,501,502,503,504
run
rendld=5.5
airyfl=t
logdfl=f
run
exit