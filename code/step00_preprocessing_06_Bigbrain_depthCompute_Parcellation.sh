#########=================depth===================###########
mripy_create_hd_mesh.ipy -f pial smoothwm  -i ../whole_SUMA -o ../whole_SUMAhd -p 3 -j 6
mripy_compute_depth_fhguoEdited.ipy -b bigbrain_SurfVol.nii  -s /Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/whole_SUMA  -l L14


#########=================HCP-MMP-1.0===================###########
3dSurf2Vol    \
-spec bigbrain_lh.spec \
-surf_A lh.smoothwm.gii \
-surf_B lh.pial.gii \
-sv bigbrain_SurfVol.nii \
-grid_parent bigbrain_1mm.nii \
-map_func max \
-f_index points  \
-f_steps 20      \
-f_p1_fr 0 -f_pn_fr 0  \
-sdata  lh.HCP-MMP-1.0.label.gii \
-prefix lh.HCP-MMP-1.0.nii  \
-overwrite

3dSurf2Vol    \
-spec bigbrain_rh.spec \
-surf_A rh.smoothwm.gii \
-surf_B rh.pial.gii \
-sv bigbrain_SurfVol.nii \
-grid_parent bigbrain_1mm.nii \
-map_func max \
-f_index points  \
-f_steps 20      \
-f_p1_fr 0 -f_pn_fr 0   \
-sdata  rh.HCP-MMP-1.0.label.gii \
-prefix rh.HCP-MMP-1.0.nii  \
-overwrite

3dcalc -prefix Atlas_HCP-MMP-1.0.nii -a rh.HCP-MMP-1.0.nii -b lh.HCP-MMP-1.0.nii  -expr '(a+180)*step(a)*step(-b+0.1) + (b-180)*step(b)' -overwrite

freeview -f lh.pial.gii:annot=lh.HCP-MMP-1.0.annot



## load laminar profile
3dresample  -master bigbrain_SurfVol.nii   -rmode NN  -prefix rm.atlas.nii   -input Atlas_HCP-MMP-1.0.nii   -overwrite
outputDir='/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/whole_SUMA/LaminarProfile_HCP_MMP1p0'
numlayers=6
step=1 / $numlayers
step=$(echo "scale=5; 1 / $numlayers" | bc)
for i in $(seq 1 $numlayers); do
    upb=$(echo "scale=5; $step * $i + 0.0001" | bc)
    downb=$(echo "scale=5; $step * ($i - 1)" | bc)
    echo $upb
    echo $downb

    3dcalc -prefix rm.depth${i}.nii -a bigbrain_SurfVol.depth.nii -expr 'step(a-'$downb')*step('$upb'-a)' -overwrite
done


NumROIs=360
for i in $(seq 1 $NumROIs); do
    uproi=$(echo "scale=5; $i + 0.5" | bc)
    downroi=$(echo "scale=5; $i - 0.5" | bc)
    3dcalc -prefix rm.ROI.nii -a rm.atlas.nii -expr 'step(a-'$downroi')*step('$uproi'-a)' -overwrite
    for j in $(seq 1 $numlayers); do
        3dcalc -prefix rm.tmpROI.nii -a rm.ROI.nii -b rm.depth${j}.nii -expr 'step(a)*step(b)' -overwrite
        3dmaskave -q -mask rm.tmpROI.nii bigbrain_SurfVol.nii > ${outputDir}/ROI${i}_depth${j}.1D
    done
done
rm -f rm*


## load V1 laminar profile
3dresample  -master bigbrain_SurfVol.nii   -rmode NN  -prefix rm.atlas.nii   -input Atlas_HCP-MMP-1.0.nii   -overwrite
outputDir='/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/whole_SUMA/LaminarProfile_HCP_MMP1p0_V1'
numlayers=50
step=$(echo "scale=5; 1 / $numlayers" | bc)

3dcalc -prefix rm.ROI.nii -a rm.atlas.nii -expr 'step(a-0.5)*step(1.5-a)' -overwrite
for i in $(seq 1 $numlayers); do
    upb=$(echo "scale=5; $step * $i + 0.0001" | bc)
    downb=$(echo "scale=5; $step * ($i - 1)" | bc)
    echo $upb
    echo $downb

    3dcalc -prefix rm.tmpROI.nii -a bigbrain_SurfVol.depth.nii -b rm.ROI.nii -expr 'step(b)*step(a-'$downb')*step('$upb'-a)' -overwrite
    3dmaskave -q -mask rm.tmpROI.nii bigbrain_SurfVol.nii > ${outputDir}/V1_depth${i}.1D
    rm -f rm.tmpROI.nii
done

rm -f rm*





## load laminar profile(on server)
3dresample  -master bigbrain_SurfVol.nii   -rmode NN  -prefix rm.atlas.nii   -input Atlas_HCP-MMP-1.0.nii   -overwrite
outputDir='/ifs/loni/groups/loft/FanhuaGuo/Data/BigBrain/LaminarProfile_HCP_MMP1_layer50'
numlayers=50
step=1 / $numlayers
step=$(echo "scale=5; 1 / $numlayers" | bc)
for i in $(seq 1 $numlayers); do
    upb=$(echo "scale=5; $step * $i + 0.0001" | bc)
    downb=$(echo "scale=5; $step * ($i - 1)" | bc)
    echo $upb
    echo $downb

    3dcalc -prefix rm.depth${i}.nii -a bigbrain_SurfVol.depth.nii -expr 'step(a-'$downb')*step('$upb'-a)' -overwrite
done


NumROIs=360
for i in $(seq 1 $NumROIs); do
    uproi=$(echo "scale=5; $i + 0.5" | bc)
    downroi=$(echo "scale=5; $i - 0.5" | bc)
    3dcalc -prefix rm.ROI.nii -a rm.atlas.nii -expr 'step(a-'$downroi')*step('$uproi'-a)' -overwrite
    for j in $(seq 1 $numlayers); do
        3dcalc -prefix rm.tmpROI.nii -a rm.ROI.nii -b rm.depth${j}.nii -expr 'step(a)*step(b)' -overwrite
        3dmaskave -q -mask rm.tmpROI.nii bigbrain_SurfVol.nii > ${outputDir}/ROI${i}_depth${j}.1D
    done
done
rm -f rm*


