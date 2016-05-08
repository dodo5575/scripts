#display update off

source /home/cli56/scripts/viewProc.tcl

set psf windmillfull2.psf
set pdb windmillfull2.pdb
#set dcd3 /data/shared/cli56/square2plate_ASiO2_LH_1MKCl_XnotPeriodic/square2plate_ASiO2_LH_1MKCl_XnotPeriodic_500mV/output/square2plate_ASiO2_LH_1MKCl_XnotPeriodic_500mV.1_21_stride_100_Wrap.dcd
set id3 [mol load psf $psf pdb $pdb]
reSet
display resetview

color Resname ADE blue2
color Resname THY green2
color Resname CYT red2
color Resname GUA orange3


mol addrep $id3
mol modcolor 0 $id3 SegName
mol modstyle 0 $id3 NewCartoon 0.900000 10.000000 1.000000 0
#mol modstyle 0 $id3 Surf 1.00000 0.000000
mol modselect 0 $id3 backbone
mol modmaterial 0 $id3 AOChalky
#mol showperiodic $id3 0 yY
#mol numperiodic $id3 0 6
#mol drawframes $id3 0 {0}

mol addrep $id3
mol modcolor 1 $id3 SegName
mol modstyle 1 $id3 Licorice 0.300000 50.000000 50.000000
mol modselect 1 $id3 not backbone and noh
mol modmaterial 1 $id3 AOChalky
#mol showperiodic $id3 1 yY
#mol numperiodic $id3 1 6
#mol smoothrep $id3 1 5

#mol addrep $id3
#mol modcolor 2 $id3 ColorID 19
#mol modstyle 2 $id3 QuickSurf 1.950000 7.609375 3.900000 1.000000
#mol modselect 2 $id3 resname DFPC and y < 0
#mol modmaterial 2 $id3 Transparent
#mol showperiodic $id3 2 yY
#mol numperiodic $id3 2 6
#mol smoothrep $id3 2 5

#mol addrep $id3
#mol modcolor 3 $id3 ColorID 0
#mol modstyle 3 $id3 QuickSurf 2.610000 17.996582 2.610000 1.000000
#mol modselect 3 $id3 water and y < 0
#mol modmaterial 3 $id3 water
##mol showperiodic $id3 3 yY
##mol numperiodic $id3 3 6
##mol smoothrep $id3 3 5
#
#mol addrep $id3
#mol modcolor 4 $id3 ColorID 1
#mol modstyle 4 $id3 Licorice 0.800000 50.000000 50.000000
#mol modselect 4 $id3 resname 3CHL 5CHL and noh
#mol modmaterial 4 $id3 Opaque
##mol showperiodic $id3 4 yY
##mol numperiodic $id3 4 6
##mol smoothrep $id3 4 5

#mol addrep $id3
#mol modcolor 5 $id3 ColorID 0
#mol modstyle 5 $id3 QuickSurf 2.500000 0.500000 1.000000 1.000000
#mol modselect 5 $id3 water and y < 0
#mol modmaterial 5 $id3 water
##mol showperiodic $id3 3 yY
##mol numperiodic $id3 3 6

set all [atomselect top all]
#$all frame 0

#animate goto 0
display resize 800 700
#display resize 600 450
display resetview
rotate y by -45
rotate x by 45
scale by 1.65
translate by 0.05 0.05 0
#rotate y by -5
#rotate x by 15
#rotate z by 1

#display projection Perspective
display antialias on
display shadows on
display ambientocclusion on
display aoambient 0.6
display aodirect 0.8
display depthcue on
display cuemode linear
display cuestart 1.500000
display cueend 2.7500000
#display backgroundgradient on
#color Display BackgroundBot red2
#color Display BackgroundTop blue2

#render aasamples TachyonInternal 32
#render aosamples TachyonInternal 32
render aasamples TachyonLOptiXInternal 100 
render aosamples TachyonLOptiXInternal 100

#set nFrames [molinfo top get numframes]

#animate goto 235
#display update on

#display dof on
#display dof_fnumber 47.0  
#display dof_focaldist 0.84

#render TachyonInternal funnel_WI_y.tga
render TachyonLOptiXInternal windmillfull2.tga

exit
