input=ch2xhe_para.inp
job=ch2x_p1
ener=913.693
jtot2=160
fstfac=1.10
tolai=1.02
spac=0.03
jmax=14
show
run
printc
quit
