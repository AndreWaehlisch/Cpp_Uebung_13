set terminal png 

set style data linespoints

set output "lorentz_1.png"
splot "lorentz_1.txt" u 1:2:3

set output "lorentz_2.png"
splot "lorentz_2.txt" u 1:2:3

set output "lorentz_3.png"
splot "lorentz_3.txt" u 1:2:3
