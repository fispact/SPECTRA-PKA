

set colors classic
##lw 3 - divide lw 1 by 2
#set lt 1 dt 1 lw 0.8
#set lt 2 dt (5,2.5)
#set lt 3 dt (2.5,3)
#set lt 4 dt (1.5,2)
#set lt 5 dt (7,2,1,2)
#set lt 6 dt (3,2.5,1,2.5)
#set lt 7 dt (2.5,2.5,2.5,7.5)
#set lt 8 dt (1,2,6,2,1,2)
#set lt 9 dt (2.5,2.5,2.5,2.5,2.5,2.5,2.5,5)

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
set logscale x
set logscale y
set logscale y2
set logscale x2
set xl "PKA energy (eV)" offset 0,1.2 font ",22"
set xr [1:3e7]
set x2r [1:3e7]

#set key off

set key nobox height 1 spacing 1 \
at graph 1.02,1.1 left top Right font ",22" \
samplen 1.3 vertical maxrows 13 
#width 0.2 samplen 1.0 horizontal maxcols 1 opaque 
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
#set border lt -1 lw 1.5 
set size square
set origin 0.0,0.3 

set yl "PKAs s^{-1} cm^{-3}" offset 2.3,0 font ",22"
set term postscript portrait color enhanced font ",25"
set output "ZR_elemental.ps"

set yr [1e6 :1.00E+12]
set y2r [1e6:1.00E+12] 
set size ratio 0.7
pl  \
"../ZR.out" index 308:308 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*6.506/91.22364157600980) \
w l lt '8' lw 5 lc rgbcolor "web-green" \
t "p",\
"../ZR.out" index 309:309 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*6.506/91.22364157600980) \
w l lt '7' lw 5 lc rgbcolor "magenta" \
t "{/Symbol a}",\
"../ZR.out" index 312:312 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*6.506/91.22364157600980) \
w l lt '1' lw 5 lc rgbcolor "black" \
t "Sr",\
"../ZR.out" index 313:313 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*6.506/91.22364157600980) \
w l lt '2' lw 5 lc rgbcolor "red" \
t "Y",\
"../ZR.out" index 314:314 \
 u (($1+$2)*1e6/2):($3*atoms_per_mole*6.506/91.22364157600980) \
w l lt '3' lw 5 lc rgbcolor "navy" \
t "{/Times-Italic Zr}"
unset output
