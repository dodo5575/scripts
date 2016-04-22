BEGIN {
  tim=0.0
}
{
  if ($1 == "END") print $1,tim
  else {
    print $1,$3,$4,$5,$6
    tim=$2
  }
}
