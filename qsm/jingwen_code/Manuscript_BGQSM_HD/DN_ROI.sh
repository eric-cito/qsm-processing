3dcalc -a Seg_ANTS.nii.gz -expr 'equals(a,17)' -prefix tmp_DNL.nii.gz
3dmask_tool -input tmp_DNL.nii.gz -prefix tmp2_DNL.nii.gz -dilate_input 2
3dcalc -a tmp2_DNL.nii.gz -b QSM_iLSQR_meanEcho_reg.nii.gz -expr 'step(b-0.03)*a' -prefix tmp3_DNL.nii.gz
3dmask_tool -input tmp3_DNL.nii.gz -prefix tmp4_DNL.nii.gz -dilate_input -1 1

3dcalc -a Seg_ANTS.nii.gz -expr 'equals(a,18)' -prefix tmp_DNR.nii.gz
3dmask_tool -input tmp_DNR.nii.gz -prefix tmp2_DNR.nii.gz -dilate_input 2
3dcalc -a tmp2_DNR.nii.gz -b QSM_iLSQR_meanEcho_reg.nii.gz -expr 'step(b-0.03)*a' -prefix tmp3_DNR.nii.gz
3dmask_tool -input tmp3_DNR.nii.gz -prefix tmp4_DNR.nii.gz -dilate_input -1 1

3dcalc -a Seg_ANTS.nii.gz -b tmp4_DNL.nii.gz -c tmp4_DNR.nii.gz -expr 'a*step(16.5-a)+17*b+18*c' -prefix Seg_ANTS_manual.nii.gz