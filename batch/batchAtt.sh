#!/bin/bash

#===============================================================================
# Grid Engine steering options
#===============================================================================
## Set the name
#$ -N computeAtt

## Submit job under TREND group
#$ -P P_trend

## Merge the stdout et stderr in a single file
#$ -j y

## Files .e et .o copied to current working directory
## -cwd

## Notify stop and kill signals before issuing them.
## -notify

## CPU time
#$ -l ct=12:00:00

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

outdir=
python /pbs/throng/trend/soft/sim/GRANDsim/simulations/computeAttenuation.py $1
cp *.att /sps/hep/trend/omartino/production/GRAND/.
