


b = 6
sample = 2

set yrange[0:0.22]
set xrange[0:1]
do for [s=1:sample] {
    set terminal gif animate delay 20 size 1240,1024
    filename = sprintf("spectrum_size18_sample%d.gif", s)
    set output filename
    do for [i=0:b] {
        
        file = sprintf("mean_spectrum_block%d_size18_sample%d.dat", i, s)
        title_file = sprintf("Mean Spectrum sample %d, for t = %d", s, 2**i * 4 * 64)
        set title title_file
        unset key
        set grid
        set xlabel "ω"
        set ylabel "I(ω)"
        set palette defined(0.4 "blue", 1.1 "red")
        plot file u 1:2:4 w l linecolor palette
    }
    unset terminal
}