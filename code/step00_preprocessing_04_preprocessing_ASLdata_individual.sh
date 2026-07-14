#########=================connect afni and suma===================###########
subj=S01
input_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/SUMA'
cd ${input_dir}/${subj}/CBF

cd ${input_dir}/${subj}
@SUMA_Make_Spec_FS -sid ${subj} 
cd ${input_dir}/${subj}/CBF
mripy_create_hd_mesh.ipy -i ../SUMA -o ../SUMAhd -p 2 -j 6



3dresample  -dxyz 0.5 0.5 0.5   -rmode Cu  -prefix rT1.nii  -input T1.nii  -overwrite
3dAutobox -prefix rT1.nii -input rT1.nii -noclust -overwrite
mripy_align_S2E.ipy -b rT1.nii  -s ../SUMAhd  -o ${subj}_SurfVol_Alnd_Exp.nii
mripy_compute_depth.ipy -b rT1.nii  -s ../SUMAhd -v ${subj}_SurfVol_Alnd_Exp.nii -l L14

###============ align control to cutrT1.nii with session 1
3dTstat -prefix CBF_v1.nii -overwrite -mean CBFv1_run1.nii
rm -f CBFv1_run1.nii

3dSkullStrip -orig_vol -prefix control1_ns.nii -overwrite -input control1.nii -overwrite
prefix=control1_al
fixed=T1.nii
moving=control1_ns.nii
base_mask=none
in_mask=none
base_mask_SyN=none
in_mask_SyN=none
antsRegistration -d 3 --float 1 --verbose \
      --output [ ${prefix}_, ${prefix}_fwd_warped.nii.gz, ${prefix}_inv_warped.nii.gz ] \
      --interpolation LanczosWindowedSinc \
      --collapse-output-transforms 1 \
      --initial-moving-transform [ ${fixed}, ${moving}, 1 ]  \
      --winsorize-image-intensities [0.005,0.995] \
      --use-histogram-matching \
      --transform translation[ 0.1 ] \
          --metric mattes[ ${fixed}, ${moving}, 1, 32, regular, 0.3 ] \
          --convergence [ 1000x300x100, 1e-6, 10 ]  \
          --smoothing-sigmas 4x2x1vox  \
          --shrink-factors 8x4x2 \
          --masks [ ${base_mask}, ${in_mask} ] \
      -t rigid[ 0.1 ] \
          -m mattes[ ${fixed}, ${moving}, 1, 32, regular, 0.3 ] \
          -c [ 1000x300x100, 1e-6, 10 ]  \
          -s 4x2x1vox  \
          -f 4x2x1  \
          -x [ ${base_mask}, ${in_mask} ] \
      -t affine[ 0.1 ] \
          -m mattes[ ${fixed}, ${moving}, 1, 32, regular, 0.3 ] \
          -c [ 1000x300x100, 1e-6, 10 ]  \
          -s 2x1x0vox  \
          -f 4x2x1  \
          -x [ ${base_mask}, ${in_mask} ]

ants2afniMatrix.py -i control1_al_0GenericAffine.mat -o antsAffine_control1.1D

3dAllineate -base rT1.nii \
			-input control1.nii \
			-interp NN  \
			-1Dmatrix_apply antsAffine_control1.1D \
			-prefix control1_al.nii -overwrite 
3dAllineate -base rT1.nii \
			-input CBF_v1.nii \
			-interp NN  \
			-1Dmatrix_apply antsAffine_control1.1D \
			-prefix rCBFv1.nii -overwrite 
3dcalc -prefix rCBFv1.nii -a rCBFv1.nii -expr 'a*step(200-a)*step(a+200)' -overwrite


###============ align control to cutrT1.nii with session 2
3dTstat -prefix CBF_v2.nii -overwrite -mean CBFv2_run1.nii
rm -f CBFv2_run1.nii

3dSkullStrip -orig_vol -prefix control2_ns.nii -overwrite -input control2.nii -overwrite
prefix=control2_al
fixed=T1.nii
moving=control2_ns.nii
base_mask=none
in_mask=none
base_mask_SyN=none
in_mask_SyN=none
antsRegistration -d 3 --float 1 --verbose \
      --output [ ${prefix}_, ${prefix}_fwd_warped.nii.gz, ${prefix}_inv_warped.nii.gz ] \
      --interpolation LanczosWindowedSinc \
      --collapse-output-transforms 1 \
      --initial-moving-transform [ ${fixed}, ${moving}, 1 ]  \
      --winsorize-image-intensities [0.005,0.995] \
      --use-histogram-matching \
      --transform translation[ 0.1 ] \
          --metric mattes[ ${fixed}, ${moving}, 1, 32, regular, 0.3 ] \
          --convergence [ 1000x300x100, 1e-6, 10 ]  \
          --smoothing-sigmas 4x2x1vox  \
          --shrink-factors 8x4x2 \
          --masks [ ${base_mask}, ${in_mask} ] \
      -t rigid[ 0.1 ] \
          -m mattes[ ${fixed}, ${moving}, 1, 32, regular, 0.3 ] \
          -c [ 1000x300x100, 1e-6, 10 ]  \
          -s 4x2x1vox  \
          -f 4x2x1  \
          -x [ ${base_mask}, ${in_mask} ] \
      -t affine[ 0.1 ] \
          -m mattes[ ${fixed}, ${moving}, 1, 32, regular, 0.3 ] \
          -c [ 1000x300x100, 1e-6, 10 ]  \
          -s 2x1x0vox  \
          -f 4x2x1  \
          -x [ ${base_mask}, ${in_mask} ]

ants2afniMatrix.py -i control2_al_0GenericAffine.mat -o antsAffine_control2.1D

3dAllineate -base rT1.nii \
			-input control2.nii \
			-interp NN  \
			-1Dmatrix_apply antsAffine_control2.1D \
			-prefix control2_al.nii -overwrite 
3dAllineate -base rT1.nii \
			-input CBF_v2.nii \
			-interp NN  \
			-1Dmatrix_apply antsAffine_control2.1D \
			-prefix rCBFv2.nii -overwrite 
3dcalc -prefix rCBFv2.nii -a rCBFv2.nii -expr 'a*step(200-a)*step(a+200)' -overwrite



###============ T1 correction
3dresample  -master rT1.nii   -rmode Cu  -prefix rT1val.nii  -input T1_val.nii  -overwrite

3dcalc -prefix rCBFcv1.nii  \
		-a rCBFv1.nii  \
		-b rT1val.nii  \
		-expr '8690.2 * exp(0.28/b) * a / (9288 * b * (1 - exp(-1.5/b)))'

3dcalc -prefix rCBFcv2.nii  \
		-a rCBFv2.nii  \
		-b rT1val.nii  \
		-expr '8690.2 * exp(0.28/b) * a / (9288 * b * (1 - exp(-1.5/b)))'





#########================= get atlas mask and change to 1D file ===================###########
##===get atlas mask from MNI152
HCPatlas_dir='/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/MNI_Glasser_HCP_2019_v1.0'
analysis_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/SUMA'
output_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data'
# subjs=('S13' 'S15' 'S16' 'S17' 'S19' 'S20' 'S21' 'S22')
subjs=('S18')
session=1
for subj in ${subjs[@]} ; do
	cd ${analysis_dir}/${subj}/CBF
	mkdir ${output_dir}/${subj}
	rm -f temp.1D

	3dNwarpApply -master rT1.nii \
				-source ${HCPatlas_dir}/HCP_atlas.nii \
				-interp NN \
				-iwarp -nwarp "T1warp_whole_1Warp.nii.gz antsAffine_whole.1D" \
				-prefix Atlas_HCP_MMP1.nii -overwrite

	3dcalc -prefix rm.mask.nii -a Atlas_HCP_MMP1.nii -expr 'step(a)' -overwrite
	3dTcat -prefix rm.out.nii -overwrite Atlas_HCP_MMP1.nii rcontrol.depth.nii rCBFcv${session}.nii
	3dmaskdump  -xyz  -o temp.1D  -mask rm.mask.nii  -noijk  rm.out.nii -overwrite
	cp  temp.1D  ${output_dir}/${subj}/Data.HCPMMP1.depth.CBF.s${session}.1D
	rm -f temp.1D
	rm -f rm*
done





##===get atlas mask from std141
# generate the atlas mask
HCPatlas_dir='/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/MNI_Glasser_HCP_2019_v1.0'
analysis_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/SUMA'
output_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data'
subjs=('S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S13' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S32' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46')
# subjs=('S18')
session=1
for subj in ${subjs[@]} ; do
	cd ${analysis_dir}/${subj}/CBF
	for hemi in lh rh ; do
        3dSurf2Vol    \
            -spec ../SUMA/std.141.${subj}_${hemi}.spec \
            -surf_A  std.141.${hemi}.white.asc \
            -surf_B  std.141.${hemi}.pial.asc \
            -sv ${subj}_SurfVol_Alnd_Exp.nii \
            -grid_parent control1.nii \
            -map_func max \
            -f_index points  \
            -f_steps 10      \
            -f_p1_fr -0 -f_pn_fr 0   \
            -sdata  ${HCPatlas_dir}/${hemi}.std.141.Glasser_HCP.lbl.niml.dset \
            -prefix rm.${hemi}.glasser.nii  \
            -overwrite
    done
    3dcalc -prefix rm.nonmix.mask.nii -a rm.lh.glasser.nii -b rm.rh.glasser.nii -expr 'step(step(a)+step(b)-2*step(a)*step(b))' -overwrite
    3dcalc -prefix rm.glasser.nii -a rm.lh.glasser.nii -b rm.rh.glasser.nii -c rm.nonmix.mask.nii  \
            -expr 'a*c + b*c + 180*step(b)*c' -overwrite
    3dROIMaker -inset rm.glasser.nii -refset rm.glasser.nii -inflate 2 -prefix rm.glasser.inf2 -overwrite
    3dTcat -prefix Atlas_HCP_MMP1_std141.nii -overwrite rm.glasser.inf2_GMI+orig.
    3dresample  -master rcontrol.depth.nii  -rmode NN  -prefix Atlas_HCP_MMP1_std141.nii  -input Atlas_HCP_MMP1_std141.nii  -overwrite
    rm -f rm*
done


# load CBF data
analysis_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/SUMA'
output_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data'
# subjs=('S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S13' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S32' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46')
subjs=('S01' 'S02' 'S03' 'S06' 'S08' 'S09' 'S13' 'S15' 'S16' 'S17' 'S19' 'S20' 'S21' 'S22')
# subjs=('S18')
session=2
for subj in ${subjs[@]} ; do
	cd ${analysis_dir}/${subj}/CBF
	mkdir ${output_dir}/${subj}
	rm -f temp.1D

	3dcalc -prefix rm.mask.nii -a Atlas_HCP_MMP1_std141.nii -expr 'step(a)' -overwrite
	3dTcat -prefix rm.out.nii -overwrite Atlas_HCP_MMP1_std141.nii rcontrol.depth.nii rCBFcv${session}.nii
	3dmaskdump  -xyz  -o temp.1D  -mask rm.mask.nii  -noijk  rm.out.nii -overwrite
	cp  temp.1D  ${output_dir}/${subj}/Data.std141HCPMMP1.depth.CBF.s${session}.1D
	rm -f temp.1D
	rm -f rm*
done


# load T1w/T1v data
analysis_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/SUMA'
output_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data'
subjs=('S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S13' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S32' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46')
subjs=('S02')
for subj in ${subjs[@]} ; do
	cd ${analysis_dir}/${subj}/CBF
	mkdir ${output_dir}/${subj}
	rm -f temp.1D

	3dcalc -prefix rm.mask.nii -a Atlas_HCP_MMP1_std141.nii -expr 'step(a)' -overwrite
    3dresample  -master rcontrol.depth.nii  -rmode NN  -prefix rm.T1w.nii  -input T1.nii  -overwrite
    3dresample  -master rcontrol.depth.nii  -rmode NN  -prefix rm.T1v.nii  -input T1_val.nii  -overwrite
	3dTcat -prefix rm.out.nii -overwrite Atlas_HCP_MMP1_std141.nii rcontrol.depth.nii rm.T1w.nii rm.T1v.nii
	3dmaskdump  -xyz  -o temp.1D  -mask rm.mask.nii  -noijk  rm.out.nii -overwrite
	cp  temp.1D  ${output_dir}/${subj}/Data.std141HCPMMP1.depth.T1.1D
	rm -f temp.1D
	rm -f rm*
done




##===delete 1D file after processing
output_dir='/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data'
subjs=('S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S13' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S32' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46')
for subj in ${subjs[@]} ; do
	cd ${output_dir}/${subj}
	# rm -f Data.HCPMMP1.depth.CBF.s1.1D
    # rm -f Data.HCPMMP1.depth.CBF.s2.1D
    # rm -f Data.std141HCPMMP1.depth.CBF.s1.1D
    # rm -f Data.std141HCPMMP1.depth.CBF.s2.1D
    rm -f Data.std141HCPMMP1.depth.T1.1D
done



