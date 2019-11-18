CPT=events.cpt

cat << EOF > $CPT
0 red 70 red
70 green 150 green
150 blue 700 blue
EOF

gmt begin ../figures/130_simulation_event_withid pdf
    gmt set FONT_ANNOT_PRIMARY 6p FORMAT_GEO_MAP ddd:mm
    gmt set MAP_FRAME_WIDTH 2p MAP_GRID_PEN_PRIMARY 0.25p,gray,2_2:1

    gmt set FONT_LABEL 6p,20 MAP_LABEL_OFFSET 4p
    gmt coast -JD135/35/30/40/7.0i -R110/160/20/50 -G244/243/239 -S167/194/223 -Bxafg -Byafg #-Lg85/11+o-0.3c/0.0c+c11+w1000k+f+u+l'scale'
    gmt meca  -Sd0.2c/0.05c -Z$CPT -M ../generated/raw_all_text.psmeca

    gmt colorbar -C$CPT -DjBR+w3c/0.3c+ml+o3.0c/0.0c -Bx+lDepth -By+lkm -L -S

    gmt plot -W1p,red << EOF
>
121.171098870338 23.5639471537689
116.241825003261 47.7205278935976
>
116.241825003261 47.7205278935976
152.758174996739 47.7205278935976
>
152.758174996739 47.7205278935976
147.828901129662 23.5639471537689
>
147.828901129662 23.5639471537689
121.171098870338 23.5639471537689
EOF

gmt end

rm $CPT