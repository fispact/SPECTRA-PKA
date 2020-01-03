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

#set label 2 point pt 7 ps (log10(1*1e6)/2) "1 MeV" at screen 0.05,0.05 offset 0.8,-0.1
set label 3 point pt 7 ps (log10(0.1*1e6)/2) "100 keV" at screen .2,0.05 offset 0.8,-0.1
set label 4 point pt 7 ps (log10(0.01*1e6)/2) "10 keV" at screen 0.35,0.05 offset 0.8,-0.1
set label 5 point pt 7 ps (log10(0.001*1e6)/2) "1 keV" at screen 0.5,0.05 offset 0.8,-0.1
set label 6 point pt 7 ps (log10(0.0001*1e6)/2) "100 eV" at screen 0.65,0.05 offset 0.8,-0.1

set border 4095
#unset border
set label 1 "(nm)" at graph -0.15,-0.1,-0.08

set key off

set term postscript portrait enhanced color 
set size 0.95

i=100




set label 20 at screen 0.6,0.85 sprintf("%.1f seconds",(i)) 

set output sprintf("pkas%is.ps",i)



file="../config_events.pka"

spl file u ($11>0?(strcol(8) eq "Zr"?\
($10<=i?$2/10:1/0):1/0):1/0):($3/10):($4/10):(log10($11*1e6)/2) \
w p   pt 7 ps variable lc rgbcolor "red", \
 




unset output
