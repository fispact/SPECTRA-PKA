

set colors classic
#lw 5 
set lt 1 dt 1 lw 0.8
set lt 2 dt (2.5,1.25)
set lt 3 dt (1.25,1.5)
set lt 4 dt (0.75,1.0)
set lt 5 dt (3.5,1,0.5,1)
set lt 6 dt (1.5,1.25,0.5,1.25)
set lt 7 dt (1.25,1.25,1.25,3.75)
set lt 8 dt (0.5,1,3,1,0.5,1)
set lt 9 dt (1.25,1.25,1.25,1.25,1.25,1.25,1.25,2.5)

atoms_per_mole=6.02214129e23
density=6.506
atomic_mass=91.22364157600980
set logscale x
set logscale y
set logscale y2
set logscale x2
set xl "PKA energy (eV)" offset 0,1.2 font ",22"
set xr [1:3e7]
set x2r [1:3e7]
#set key at graph 1.15,1.15 top left width -3 samplen 2.5 horizontal  maxcols 1


set key nobox height 1 spacing 1.0 \
at graph -0.1,-0.22 left top Right font ",18" \
width -0.5 samplen 1.0 horizontal maxcols 5 
#vertical maxrows 5  
set ytics 100 format "10^{%T}" nomirror out scale 2.0,1.0 offset 1.4,0.2 font ",20"
set xtics 100 format "10^{%T}" out nomirror scale 2.0,1.0 offset 0.0,0.6 font ",20"
set y2tics 10 format "" in mirror
set my2tics 10
set mytics 20
set x2tics in mirror tc rgbcolor "white" scale 2.0,1.0 font "Times-Roman,3" format " "
set mx2tics 10
set mxtics 20
set xtics add (10 1)
#set grid x2tics lt 2 lw 1 lc rgbcolor "grey"
#set grid ytics lt 2 lw 1 lc rgbcolor "grey"
#set border lw 2 
set size square
set origin 0.0,0.3 

set yl "PKAs s^{-1} cm^{-3}" offset 2.3,0 font ",22"
set term postscript portrait color enhanced font ",25"
set output "ZR_isotope.ps"

set yr [1e6 :1.00E+12 ]
set y2r [1e6:1.00E+12 ] 
set size ratio 0.7
pl  \
"../ZR.out" index 508:508 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '8' lw 5 lc rgbcolor "web-green" \
t "H1",\
"../ZR.out" index 509:509\
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '7' lw 5 lc rgbcolor "magenta" \
t "He4",\
"../ZR.out" index 514:514 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '1' lw 5 lc rgbcolor "black" \
t "{/Times-Italic Zr90}",\
"../ZR.out" index 513:513\
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '2' lw 5 lc rgbcolor "red" \
t "Sr88",\
"../ZR.out" index 517:517 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '3' lw 5 lc rgbcolor "navy" \
t "Sr89",\
"../ZR.out" index 519:519 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '5' lw 5 lc rgbcolor "web-blue" \
t "{/Times-Italic Zr91}",\
"../ZR.out" index 520:520 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '3' lw 5 lc rgbcolor "dark-grey" \
t "Y91",\
"../ZR.out" index 524:524 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '1' lw 5 lc rgbcolor "dark-goldenrod" \
t "{/Times-Italic Zr92}",\
"../ZR.out" index 529:529 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '6' lw 5 lc rgbcolor "cyan" \
t "Zr93",\
"../ZR.out" index 533:533 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '3' lw 5 lc rgbcolor "light-green" \
t "{/Times-Italic Zr94}",\
"../ZR.out" index 538:538 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '5' lw 5 lc rgbcolor "purple" \
t "Zr95",\
"../ZR.out" index 540:540 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '8' lw 5 lc rgbcolor "steelblue" \
t "{/Times-Italic Zr96}",\
"../ZR.out" index 543:543 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*density/atomic_mass) \
w l lt '7' lw 5 lc rgbcolor "dark-cyan" \
t "Zr97"
unset output
