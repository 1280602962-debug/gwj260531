load /mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0/analysis/eh40_23_vs_anchors_v1/poses_pdb/receptor_3RCD.pdb, receptor
hide everything, receptor
show cartoon, receptor
color gray80, receptor
load /mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0/analysis/eh40_23_vs_anchors_v1/poses_pdb/3RCD_EH40_01_mode03.pdb, EH40_01
show sticks, EH40_01
color green, EH40_01
load /mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0/analysis/eh40_23_vs_anchors_v1/poses_pdb/3RCD_EH40_02_mode01.pdb, EH40_02
show sticks, EH40_02
color cyan, EH40_02
load /mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0/analysis/eh40_23_vs_anchors_v1/poses_pdb/3RCD_EH40_23_mode02.pdb, EH40_23
show sticks, EH40_23
color magenta, EH40_23
zoom receptor
png /mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0/analysis/eh40_23_vs_anchors_v1/figures/overlay_3RCD.png, dpi=200
