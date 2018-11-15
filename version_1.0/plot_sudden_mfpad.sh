#!/bin/bash

gnf=gnufile.gp

lw1=2
lw2=2
lw3=2
lw4=2
lw5=2
lw6=2
lw7=2
lw8=2
color1=black
color2=red
color3=blue
color4=green
color5=purple
color6=gold
color7=magenta
color8=pink
lt1=1
lt2=1
lt3=1
lt4=1
lt5=1
lt6=1
lt7=1
lt8=1
dt1=1
dt2=1
dt3=1
dt4=1
dt5=1
dt6=1
dt7=1
dt8=1

generating_energy=1.65 #eV
harmonic_order=11
energy=$(awk "BEGIN {print ${generating_energy}*${harmonic_order}}" ) #eV
duration_FWMH=$(awk "BEGIN {print 20*(27.211/${energy})*0.02418884}" ) #fs
duration_sigma=$(awk "BEGIN {print ${duration_FWMH}/(2*((2*log(2))**0.5)*0.02418884) }") #au
bandwidth=$(awk "BEGIN {print 27.211/${duration_sigma}}" )

echo "Exciting pulse centered at ${generating_energy} eV. Taking the ${harmonic_order}th harmonic" 
echo "Pulse centered at ${energy} eV."
echo "Sigma = ${duration_sigma} au => Sigma energy = ${bandwidth} eV"

cat > $gnf << MAFG
set terminal x11 enhanced font "Times-Roman, 20 pts"

#set terminal postscript enhanced color font "Times-Roman, 20 pts"
#set output "cs_field_modulated_sudden_pump.eps" 

set xrange [0:30]
set style line 1 lt $lt1 lc rgb "${color1}" lw $lw1 dt $dt1 
set style line 2 lt $lt2 lc rgb "${color2}" lw $lw2 dt $dt2 
set style line 3 lt $lt3 lc rgb "${color3}" lw $lw3 dt $dt3 
set style line 4 lt $lt4 lc rgb "${color4}" lw $lw4 dt $dt4 
set style line 5 lt $lt5 lc rgb "${color5}" lw $lw5 dt $dt5 
set style line 6 lt $lt6 lc rgb "${color6}" lw $lw6 dt $dt6 
set style line 7 lt $lt7 lc rgb "${color7}" lw $lw7 dt $dt7 
set style line 8 lt $lt8 lc rgb "${color8}" lw $lw8 dt $dt8 

plot "cs_0_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+7.330961769)**2/(${bandwidth}**2))) w l ls 1 t "GS"\\
,"cs_1_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+4.209362)**2/(${bandwidth}**2))) w l ls 2 t "1S"\\
,"cs_2_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+1.99450336)**2/(${bandwidth}**2))) w l ls 3 t "2S"\\
,"cs_3_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+1.62907)**2/(${bandwidth}**2))) w lp ls 4 t "3S"\\
,"cs_4_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+1.425187)**2/(${bandwidth}**2))) w lp ls 5 t "4S"\\
,"cs_5_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+0.986554)**2/(${bandwidth}**2))) w l ls 6 t "5S"\\
,"cs_6_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+0.826926)**2/(${bandwidth}**2))) w l ls 7 t "6S"\\
,"cs_9_0.txt" u 1:(\$2)*(exp(-((\$1)-${energy}+0.333427)**2/(${bandwidth}**2))) w l ls 8 t "9S"
MAFG

gnuplot -p "$gnf"
