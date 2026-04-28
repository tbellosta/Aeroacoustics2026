# --- Terminal ---
set terminal pdfcairo size 10cm,10cm enhanced font "Times,12"
set output 'directivity.pdf'

# --- Data range ---
stats '../bin/directivity.dat' using 2 nooutput
SPL_min = STATS_min - 3
SPL_max = STATS_max + 3

# --- Polar setup ---
set polar
set angles degrees
set size square

# --- Remove rectangular frame completely ---
unset border

# --- Draw circular frame manually ---
set object 1 circle at 0,0 size graph 0.5 \
    front lw 1.5 lc rgb "black" fillstyle empty

# --- Grid ---
set grid polar 30 lw 1 dt 2 lc rgb "#888888"

unset xtics
unset ytics

set samples 720


# --- Ranges ---
set rrange [SPL_min:SPL_max]
set trange [0:360]

# --- Radial ticks ---
set rtics 10
set rtics format "%.0f"

# --- Legend (unboxed) ---
set key top left

# --- Plot ---
plot '../bin/directivity.dat' using 1:2 with linespoints \
        lw 1.5 pt 11 ps 0.8 lc rgb "black" \
        title "FWH numerical", \
     '../bin/directivity.dat' using 1:4 with lines \
        lw 2 dt 2 lc rgb "#1f77b4" \
        title "Exact"
