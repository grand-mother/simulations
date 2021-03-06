batchpath=/pbs/throng/trend/soft/sim/GRANDsim/simulations/batch
#target=/sps/hep/trend/omartino/GRAND/voltageOutput/flat/flat_freespace
target=/sps/hep/trend/omartino/GRAND/voltageOutput/hotspot-150x67km2/HS1_freespace
jsonpath=$target/jsons
outdir=/sps/hep/trend/omartino/production/GRAND
# First build list of showers already treated
cd $outdir
rm done.txt
done=`ls *.att`
for i in $done
     do     
     ID=${i::-4}
     echo $ID
     echo $ID >> done.txt
done

# Now create json files subdirs for batches
cd $jsonpath
rm -rf batch*
if [ 1 -gt 0 ]  # 
   then
   jsons=`ls *.voltage.json`
   njsons=`ls $jsons | wc -l`
   echo $njsons "json files"
  
   nb=0
   count=100
   maxjson=9
   for json in $jsons
     do
     if [ $count -gt $maxjson ]  # Now reaching max nb of json file/subdir ==> creat enew one
       then
       let nb=nb+1
       cd $jsonpath
       subdir=batch$nb
       mkdir $subdir
       cd $subdir
       count=0
     fi
   # Add to subfolder if not processed yet 
     testID=${json::-13}
     isin=`grep $testID $outdir/done.txt | wc -l`
     if [ $isin -eq 0 ]
       then
       echo "Creating sim link of" $json "in" $subdir
       ln -s ../$json ./$json
       let count=count+1
     fi
   done
fi
#


## Now launch batches
cd $batchpath
batchdirs=`ls -d $jsonpath/batch*`
for batchdir in $batchdirs
do
  qsub -P P_trend $batchpath/batchAtt.sh $batchdir
done

