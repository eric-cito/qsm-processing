


patient=PDa428_no.consent.yet-addpost
saveTo=~/headers_$patient.json

rm $saveTo
echo "[" >> $saveTo

for loc in /data/morrison/data/parkinsons/retro_clin/$patient/qsm/qsm_dicoms/1.*/*.1.dcm;
do
	echo $loc
	dcm2json $loc | tee --append $saveTo
	echo "," >> $saveTo
done

echo "{} ]" >> $saveTo
