#!/bin/bash 
source ~/.bashrc 
echo date1: `date`
wrkdir=/workfs/ybj/jlzhang/grand
#cd $wrkdir

pool=sim/flat-100x100km2
basedir=/scratchfs/ybj/jlzhang/grand
myfields=$basedir/$pool/fields
#dirs=$(ls -d /sps/hep/trend/grand/massProd/exampleShower/??????)
myvoltages=$basedir/$pool/voltages
myseeds=$basedir/$pool/seeds

todofile=$basedir/todofile.txt

if [ ! -f $todofile ]; then
   find $myfields/ -name "*.tgz" |xargs ls > $todofile
fi
ntotal=`cat $todofile |wc -l`


#donefile=${0}.done
donefile=$basedir/computevoltage_Xmax.sh.done
echo $donefile

if [ -f $donefile ] ; then
   lastfieldid=$(tail -1 $donefile)
   echo "Last field processed: $lastfieldid "
   #ls -lhtr  $myfields/$lastfieldid/out.tar.gz

   #dirs=`find $myfields/ -name "*.tgz" -type f  -newer ${lastfieldid}.tgz |xargs ls -tr` #jlzhang, when rerun, avoid the time order bug   
   #ndirs=`find $myfields/ -name "*.tgz" -type f  -newer ${lastfieldid}.tgz | wc -l`
   lastfieldidn=`grep -n $lastfieldid $todofile  | gawk -F ':' '{print $1}'` #if all field files are ready
   dirs=`tail -$lastfieldidn $todofile` #or tail -$lastfieldidn $todofile > /tmp/${0}.txt
   ndirs=$((ntotal-n))
else
   dirs=`cat $todofile`
   ndirs=$ntotal
   #dirs=/sps/trend/jlzhang/production/grandproto/fields/10850122/out.tar.gz #for test only
fi
echo "Nb of fields to be processed: $ndirs" 
if [ $ndirs -eq 0 ] ;then
   echo "All are done!!!"
   exit
fi
echo "Press a key to proceed..."
#stty -icanon min 0 time 100 
#read touche
#stty icanon

nout=0
ntracet=0
#for dir in  $myfields/E*/*
for dir in $dirs 
#for dir in /scratchfs/ybj/jlzhang/grand/sim/flat-100x100km2/fields/E.5e18_X.181841_Y.39925_Z.541_T.91_P.33_D.12135181768782943.tgz
do

  echo "Now processing  folder " $dir
  #cd $dir
  #rundir=${dir%/*}
  rundir=${dir%.tgz}

  runname=${rundir##*/}

  ntracet=`tar ztf $dir |wc -l`
  ntracet=$((ntracet-2))
  if [ -d $rundir ] ; then  
     ntrace=`ls $rundir|wc -l`
     ntrace=$((ntrace-2))
  fi  
  #nout=`ls $myvoltages/$runname|wc -l`
  #for satety, we can use the biggst antenna id.
  #FILETEST="$myvoltages/$runname/doing.txt"
  FILETEST="$myvoltages/$runname"
  #echo $FILETEST
  #for safety, we can check whether the no. of trace and out files are equal
  #if [ $ntracet = $nout -a -f  "$FILETEST" -a ! -s "$FILETEST" ] ; then
  #if [ $ntracet = $nout -a -s $myvoltages/$runname.tgz ] ; then
  #if [ -f $FILETEST ] ; then
  #this if loop and next, this loop is true, but there may be no voltages files!!! this if loop avoid multi process doing the same run. so for the check script, comment the continue in this loop.
  if [ -d $FILETEST ] ; then
     nout=`ls $myvoltages/$runname|wc -l`
     echo "Field is being processed: $myvoltages/$runname exists"
     continue
  elif [ -f $FILETEST.tgz ]; then
       nout=`tar ztf $$FILETEST.tgz |wc -l`
       nout=$((nout-2))
       echo "Field is being processed: $myvoltages/$runname.tgz exists"
       continue
  fi

  if [ $ntracet = $nout ] ; then #note the blank space before and after "="  
     echo "Field was already processed:  $ntracet==$nout, file $myvoltages/$runname exists!"
     echo $rundir >>  $donefile
    #echo  "Writing fieldID to donefile."
    #echo $fieldid  >> $wrkdir/donefile
     continue
  fi

  mkdir -p $myvoltages/$runname
  #touch $FILETEST

  #if [ $ntracet!=$ntrace -a ! -f $myfields/$runname/a0.trace ] ; then
  if [ $ntracet!=$ntrace  ] ; then
     echo "extract $myfields/${runname}.tgz"
     cd $myfields
     rm -rf $myfields/${runname}
     tar zxf $myfields/${runname}.tgz  
  fi

  echo "Now computing antenna response for $rundir..."
  thetat=${runname#*T.}
  theta=${thetat/\_P*}
  phit=${runname#*P.}
  phi=${phit/\_D*}

  cd $myseeds
  showid_json=`grep $runname events*.json  |gawk -F "[:]" '{print $1}'`
  cd $wrkdir  #where npy in
  #echo  $wrkdir/computevoltage_Xmax.py $rundir  $myvoltages/$runname  1 1 $myseeds/$showid_json 
  python  $wrkdir/computevoltage_Xmax.py $rundir  $myvoltages/$runname  1 1 $myseeds/$showid_json 
  #python ./computevoltage_Xmax.py /scratchfs/ybj/jlzhang/grand/sim/flat-100x100km2/fields/E.9e18_X.9970_Y.8269_Z.1585_T.90_P.108_D.4959329178210273  /scratchfs/ybj/jlzhang/grand/sim/flat-100x100km2/voltages/E.9e18_X.9970_Y.8269_Z.1585_T.90_P.108_D.4959329178210273 15 1  /scratchfs/jlzhang/grand/sim/flat-100x100km2/seeds/events.29939407.11.3.json 
  pystat=$?
  if [  $pystat = 0 ] ; then
     #echo
     echo $rundir >>  $donefile
     
     rm -rf  $myfields/${runname}

     cd $myvoltages
     rm -rf $myvoltages/$runname.tgz
     tar -czf $myvoltages/$runname.tgz $runname 
  else
     echo "$wrkdir/computevoltage_Xmax.py error!!!"
     rm -rf $myfields/$runname $myvoltages/$runname $myvoltages/$runname.tgz
     continue 
  fi

done

echo "Done!"
echo date2: `date`
