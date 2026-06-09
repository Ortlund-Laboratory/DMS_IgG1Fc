load data/Fc_FcgR3a_mapping.pse

select blue_sites, resi 240+242+265+271+282+301+309+315+317+318+319+326+344+359+378+380+390+394+396+423+424+427+433+434+436+438+439
color blue, blue_sites

select red_sites, resi 233+236+237+238+239+241+251+257+267+269+270+277+293+297+298+304+310+323+325+327+329+367+374+381+392+414+429+430
color red, red_sites

show surface, surface_A
show surface, surface_B
