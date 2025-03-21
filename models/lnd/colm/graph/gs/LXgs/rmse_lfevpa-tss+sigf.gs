'reinit';
'open ../ctl/rmse_lfevpa-tss.ctl';
'open ../ctl/rmse_lfevpa-sigf.ctl';
'set t 1 204';
'set vpage 1 10 0 4.5';
'set grads off';
'set xlab off';
'set ccolor 1';
'set cstyle 1';
'set cmark  0';
'd rmse_la';
'set ylpos 0 r';
'set ccolor 1';
'set cstyle 1';
'set cmark  2';
'd ts';
'run cbar_line -x 2.1 -y 4.5 -m "0" 2 -l 1 1 -t "SD of canopy latent heat flux (W/m`a2`n)" "temperature of first soil layer (K)"'; 
'set vpage 1 10 4 8.5';
'set grads off';
'set xlab off';
'set ylpos 0 r';
'set ccolor 1';
'set cstyle 1';
'set cmark  2';
'd sf.2';
'set ylpos 0 l';
'set ccolor 1';
'set cstyle 1';
'set cmark  0';
'd rmse_la.2';
'run cbar_line -x 2.1 -y 4.5 -m "0" 2 -l 1 1 -t "SD of canopy latent heat flux (W/m`a2`n)" "fraction of veg cover"';
'enable print /g4/lix/rmse_lfevpa-tss+sigf.gmf';
'print';
'disable print';
'disable fwrite';
;
