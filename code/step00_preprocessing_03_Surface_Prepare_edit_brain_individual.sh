freeview \
	-v mri/T1.mgz \
	   mri/brainmask.mgz \
	   mri/wm.mgz \
	-f surf/lh.white:edgecolor=yellow \
	   surf/lh.pial:edgecolor=red \
	   surf/rh.white:edgecolor=yellow \
	   surf/rh.pial:edgecolor=red



###############
% just modified pial
recon-all -autorecon-pial -subjid S28_edited \
	  -sd /ifs/loni/groups/loft/FanhuaGuo/Experiment/Pulsatility_VASO2023/SUMA/Modified_Edited \
	  -parallel -openmp 12 \
	  -hires
	 
