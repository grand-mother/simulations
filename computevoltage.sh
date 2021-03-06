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

#todofile=$basedir/diff.txt #run diff.sh first
todofile=$basedir/todofile.txt
#same as run gen_todofile.sh first
if [ ! -f $todofile ]; then
   find $myfields/ -name "*.tgz" |xargs ls > $todofile
fi
ntotal=`cat $todofile |wc -l`

#donefile=${0}.done
donefile=$basedir/computevoltage.sh.done
echo $donefile


if [ -f $donefile ] ; then
   lastfieldid=$(tail -1 $donefile)
   echo "Last field processed: $lastfieldid "

   lastfieldidn=`grep -n $lastfieldid $todofile  | gawk -F ':' '{print $1}'` #if all field files are ready. if cann't grep, lastfieldidn=0
   ndirs=$((ntotal-lastfieldidn))
   dirs=`tail -$lastfieldidn $todofile ` #or tail -$lastfieldidn $todofile > /tmp/${0}.txt
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

#for dir in  $myfields/E*/*
for dir in $dirs 
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
  FILETEST="$myvoltages/$runname"
  #echo $FILETEST
  #for safety, we can check whether the no. of trace and out files are equal
  #if [ $ntracet = $nout -a -f  "$FILETEST" -a ! -s "$FILETEST" ] ; then
  #if [ $ntracet = $nout -a -s $myvoltages/$runname.tgz ] ; then
  #this if loop and next, this loop is true, but there may be no voltages files!!! this if loop avoid multi process doing the same run. so for the check script, comment this loop.
  if [ -d $FILETEST ] ; then
     nout=`ls $myvoltages/$runname|wc -l`
     echo "Field is being processed: $myvoltages/$runname exists"
     continue
  elif [ -f $FILETEST.tgz ]; then
       nout=`tar ztf $FILETEST.tgz |wc -l`
       nout=$((nout-2))
       echo "Field is being processed: $myvoltages/$runname.tgz exists"
       continue
  fi
  
  if [ $ntracet = $nout ] ; then #note the blank space before and after "="  
     echo "Field was already processed:  $ntracet==$nout, file $myvoltages/$runname exists!"
     echo $rundir >>  $donefile
    #echo  "Writing fieldID to donefile."
    continue
  fi

  mkdir -p $myvoltages/$runname
  #touch $FILETEST

  #if [ $ntracet!=$ntrace -a ! -f $myfields/$run/a0.trace ] ; then
  if [ $ntracet!=$ntrace  ] ; then
     echo "extract $myfields/${runname}.tgz"
     cd $myfields
     rm -rf $myfields/${runname}
     tar zxf $myfields/${runname}.tgz  
  fi

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
  pystat=$?
  if [ $pystat = 0  ] ; then
     #echo
     echo $rundir >>  $donefile
     
     rm -rf  $myfields/${runname}

     cd $myvoltages
     rm -rf $myvoltages/$runname.tgz
     tar -czf $myvoltages/$runname.tgz $runname 
     #rm -f $FILETEST
  else
     echo "$wrkdir/computevoltage.py error!!!"
     rm -rf $myfields/$runname $myvoltages/$runname $myvoltages/$runname.tgz
     continue
  fi

done

echo "Done!"
echo date2: `date`
