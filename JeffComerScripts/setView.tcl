set view {{{1 0 0 0.0212044} {0 1 0 7.85836} {0 0 1 1.37973} {0 0 0 1}} {{0.998918 -0.0447759 0.012562 0} {0.00678753 0.407607 0.913133 0} {-0.0460067 -0.91206 0.40747 0} {0 0 0 1}} {{0.0171213 0 0 0} {0 0.0171213 0 0} {0 0 0.0171213 0} {0 0 0 1}} {{1 0 0 -0.01} {0 1 0 -0.26} {0 0 1 -0.03} {0 0 0 1}}}
set molList [molinfo list]
foreach m $molList {
	molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $view
}



