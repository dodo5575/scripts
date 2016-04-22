Scripts written by Mahika Dubey, May-July 2013
-----------------------------------------------

setting.tcl
- This script moves the display to the initial setting for most of my videos
- First hit the '=' key to move the display to the 'home' position of the system, then source 'setting.tcl' from the Tcl/Tk console

setting2.tcl
- Another script for another initial display setting
- First hit the '=' key to change the display to the 'home' position of the system, then source 'setting2.tcl' from the Tcl/Tk console

rotatex.tcl
- This script can be used to create a VMD animation in which the display rotates around the DNA Origami
- First fix the display using 'setting.tcl'
- To run this script, source 'rotatex.tcl' from the Tcl/Tk console, then type 'enablemoviecallback.' In the VMD Movie Maker menu (VMD Main -> Extensions -> Visualization -> Movie Maker), enter all specifications and click 'Make Movie'

movie1.tcl
- This script creates a VMD animation in which the camera moves around and then inside the DNA origami bundle and rotates the view
- The initial setting is fixed by 'setting.tcl'
- To run this script, source 'movie1.tcl' from the Tcl/Tk console, then type 'enablemoviecallback.' In the VMD Movie Maker menu, enter all specifications and clock 'Make Movie'

movie2.tcl
- This script creates a VMD animation in which the camera moves into the DNA Origami bundle and rotates around inside
- The initial setting is fixed by 'setting2.tcl'
- to run this script, source 'movie2.tcl' from the Tcl/Tk console, then type 'enablemoviecallback.' In the VMD Movie Maker menu, enter all specifications and click 'Make Movie'

frame_maker.tcl
- This script generates frames for a VMD animation rendered using snapshot
- The intial setting for this script is fixed by 'setting.tcl'
- Once the display is correct, source 'frame_maker.tcl' from the Tcl/Tk console to create the frames

