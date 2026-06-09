# Load your session
load data/Fc_FcgR3a_mapping.pse

# Select surfaces
select surfA, surface_A
select surfB, surface_B

# Example: BLUE sites (replace manually or script generate)
select blue_sites, resi 220+240+250
color blue, blue_sites

# RED sites
select red_sites, resi 230+260+270
color red, red_sites

# Apply to both surfaces
show surface, surfA
show surface, surfB
