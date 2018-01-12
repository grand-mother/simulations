#!/bin/bash 
source ~/.bashrc 
echo date1: `date`
wrkdir=/workfs/ybj/jlzhang/grand/
#cd $wrkdir

pool=sim/flat-100x100km2
basedir=/scratchfs/ybj/jlzhang/grand
myfields=$basedir/$pool/fields
#dirs=$(ls -d /sps/hep/trend/grand/massProd/exampleShower/??????)
myvoltages=$basedir/$pool/voltages

#donefile=${0}.done
donefile=$basedir/computevoltage.sh.done
echo $donefile

if [ -f $donefile ] ; then
   lastfieldid=$(tail -1 $donefile)
   echo "Last field processed: $lastfieldid "
   #ls -lhtr  $myfields/$lastfieldid/out.tar.gz

   dirs=`find $myfields/ -name "*.tgz" -type f  -newer ${lastfieldid}.tgz |xargs ls -tr` #jlzhang, when rerun, avoid the time order bug   
   ndirs=`find $myfields/ -name "*.tgz" -type f  -newer ${lastfieldid}.tgz | wc -l`
else
   dirs=`find $myfields/ -name "*.tgz" |xargs ls -tr`
   ndirs=`find $myfields/ -name "*.tgz" |wc -l`
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

#for dir in  $myfields/E*/*
for dir in $dirs 
do

  echo "Now processing  folder " $dir
  #cd $dir
  #rundir=${dir%/*}
  rundir=${dir%.tgz}

  run=${rundir##*/}

  ntracet=`tar ztf $dir |wc -l`
  ntracet=$((ntracet-2))
  ntrace=`ls $rundir|wc -l`
  ntrace=$((ntrace-2))
  nout=`ls $myvoltages/$run|wc -l`
  #for satety, we can use the biggst antenna id.
  FILETEST="$myvoltages/$run/out_0.txt"
  #echo $FILETEST
  #for safety, we can check whether the no. of trace and out files are equal
  #if [ $ntracet==$nout -a -f  "$FILETEST" -a ! -s "$FILETEST" ] ; then
  if [ $ntracet==$nout -a  -s "$FILETEST" ] ; then
    echo "Field was already processed:  $ntracet==$nout, file $myvoltages/$run/out_0.txt exists!"
    #echo  "Writing fieldID to donefile."
    #echo $fieldid  >> $wrkdir/donefile
    continue
  else
    echo "No file $FILETEST"
  fi

  #mkdir -p 
  #if [ $ntracet!=$ntrace -a ! -f $myfields/$run/a0.trace ] ; then
  if [ $ntracet!=$ntrace  ] ; then
     echo "extract $myfields/${run}.tgz"
     cd $myfields
     tar zxf $myfields/${run}.tgz  
  fi
  mkdir -p $myvoltages/$run
  touch $FILETEST

  echo "Now computing antenna response for $rundir..."
  runname=${rundir##*/}
  thetat=${runname#*T.}
  theta=${thetat/\_P*}
  phit=${runname#*P.}
  phi=${phit/\_D*}

  echo theta: $theta $phi $rundir
  cd $wrkdir  #where npy in
  python  $wrkdir/computevoltage.py $theta $phi 1 $rundir  $myvoltages/$run 1
  #python computevoltage.py 70 40 1 /sps/hep/trend/grand/massProd/exampleShower/shower1 1
  #python computevoltage.py 90 32 1 /scratchfs/ybj/jlzhang/grand/showers/E.1e17/E.1e17_X.10083_Y.71097_Z.217_T.90_P.32_D.8178857730898142 1 
  if [  !$? ] ; then
     #echo
     echo $rundir >>  $donefile
     
     cd $myvoltages
     tar -czf $myvoltages/$run.tgz $run 
  fi

done

echo "Done!"
echo date2: `date`
