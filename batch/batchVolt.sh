#!/bin/bash

#===============================================================================
# Grid Engine steering options
#===============================================================================
## Set the name
#$ -N computeVolt

## Submit job under TREND group
#$ -P P_trend

## Merge the stdout et stderr in a single file
#$ -j y

## Files .e et .o copied to current working directory
## -cwd

## Notify stop and kill signals before issuing them.
## -notify

## CPU time
#$ -l ct=01:00:00

## Memory
#$ -l vmem=16.0G

## Disk space
#$ -l fsize=30.0G

## Request /sps/:
#$ -l sps=1

#$ -l irods=1

## Files .e et .o copied to the directory
#$ -e /sps/hep/trend/omartino/production/GRAND/
#$ -o /sps/hep/trend/omartino/production/GRAND/

#===============================================================================

# First fetch 3 input files
tardir=/sps/hep/trend/omartino/GRAND/voltageOutput/hotspot-150x67km2/HS1_freespace
#tardir=/sps/hep/trend/omartino/GRAND/voltageOutput/hotspot-150x67km2/HS1_ground
#attdir=/sps/hep/trend/omartino/GRAND/attFiles/flat/
attdir=/sps/hep/trend/omartino/GRAND/attFiles/HS1/
jsondir=$1  # Folder containing json files

cd $jsondir
jsons=`ls *.voltage.json`
for json in $jsons
  do
  E=${json:2:4}
  th=${json:9:2}
  if [ "${json:17:1}" = "_" ]; then # 3 digits for azimuth
    phi=${json:14:3}
    lat=${json:21:2}
    long=${json:27:2}
    ID=${json:42:5}
  elif [ "${json:16:1}" = "_" ]; then # 2 digits for azimuth
    phi=${json:14:2}
    lat=${json:20:2}
    long=${json:26:2} 
    ID=${json:41:5}
  elif [ "${json:15:1}" = "_" ]; then # 1 digit for azimuth
    phi=${json:14:1}
    lat=${json:19:2}
    long=${json:25:2} 
    ID=${json:40:5}
  fi
  cd $tardir/La$lat'_Lo'$long/E$E/Z$th/A$phi
  tarfile=`ls *$ID*.tgz`
  truetarfile=$tardir/La$lat'_Lo'$long/E$E/Z$th/A$phi/$tarfile
  ln -s $truetarfile $jsondir/$tarfile
done

cd $jsondir
# Now run job
python /pbs/throng/trend/soft/sim/GRANDsim/simulations/computeFiltResponse.py $jsondir $attdir $jsondir/

# Now copy output to result dir
cp *.json /sps/hep/trend/omartino/production/GRAND/.
