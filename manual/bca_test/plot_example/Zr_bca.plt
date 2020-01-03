box_units=500
latt=3.232

xmax=box_units*latt/10
ymax=box_units*latt*cos(pi/6)/10
zmax=box_units*latt*sqrt(3)/10

set xr [0:xmax]
set yr [0:ymax]
set zr [0:zmax]

#set key off
#set xlabel "x"
set xtics 40 out offset 0,-1
set ytics 40 out offset 1,-1
set ztics out
set xyplane at 0
set view 75,30
set view equal xyz

set border 4095
#unset border
set label 1 "(nm)" at graph -0.15,-0.1,-0.08

set key off

set term postscript portrait enhanced color 
set size 0.95

i=100




set label 2 at screen 0.6,0.85 sprintf("%.1f seconds",(i)) 

set output sprintf("bca%is.ps",i)



file="../config_bca.dat"

  
spl file u ($2<=i?$3/10:1/0):($6==1?$4/10:1/0):($5/10) \
w p lt 1 pt 7 ps 0.5 lc rgbcolor "green" t "pka interactions", \
file u ($2<=i?$3/10:1/0):($6==2?$4/10:1/0):($7<=1.0?$5/10:1/0) \
w p lt 1 pt 7 ps 0.3 lc rgbcolor "#56b4e9" t "recoils", \
file u ($2<=i?$3/10:1/0):($6==2?$4/10:1/0):($7>1.0?$5/10:1/0) \
w p lt 3 pt 7 ps 0.2 lc rgbcolor "red" t "secondary recoils"



unset output
