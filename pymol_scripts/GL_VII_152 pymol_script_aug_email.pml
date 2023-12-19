
bg_color white

load 290phage.cif.gz
show sphere
colour grey
zoom center, 500  
rotate y, 90






select glycans, chain S
set sphere_scale, 30, glycans
set sphere_transparency, 0.2000, glycans
colour orange, glycans  # was green
select junk, resi 9999  # unselect everything

ray 1600, 1200
png orange_1000.png  





