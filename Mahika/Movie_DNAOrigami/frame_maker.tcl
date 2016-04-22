set numFrames [molinfo top get numframes]

for {set frame 0} {$frame<$numFrames} {incr frame 1} {
  animate goto $frame
  if {$frame <= 71} {
      rotate y by 1.2;
      translate by 0 0 0.01;
  } elseif {$frame <= 160} {
      scale by 1.02;
  } else {
      rotate y by 1.3;
  }

  render snapshot "origami[format %05d $f].tga"
}
