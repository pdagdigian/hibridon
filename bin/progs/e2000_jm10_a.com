inp=Hch4.inp
iop=1
job=ajm10
jtot2=105
jmax=10
ener=2000
show
run
printc
trnprt
quit