inp=arn2_100.inp
noprin=f
airyfl=f
spac=0.02
rendld=25
jmax=24
jtot1=24
jtot2=24
energy=1000
job=Bench
run
spac=0.05
tolai=1.05
logdfl=f
airyfl=t
run
inp=arn2_dxsec.inp
jmax=8
jout,6,0,2,4,6,8,10
jtot2=150
run
difcrs,,8,0,8,0,0,40,.05,0,150
tenxsc,,,4,4,0,0,100,4,4
exit
