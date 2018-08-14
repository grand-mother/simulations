
#==============================================================================

# Script created to analyze Xmax GRAND simulations

#==============================================================================

SCRIPTDIR=/Users/guepin/Documents/GRAND/UHECR_Xmax/ # Folder where your scripts are saved


# GEANTDIR=/home/zilles/Geant4-SKA/LORA

# How to proceed with numbering the simulations is not consistent so far. If there is one number missing the script just skips it.
# For running the analysis script you have nevertheless to hand over a list file with the simuations you like to use

#NR=CR170_83deg
#filename=SimlistTest170_83.txt

#NR=CR175_83deg
#filename=SimlistTest175_83.txt


#NR=CR180_83deg
#filename=SimlistTest180_83.txt

#NR=CR185_83deg
#filename=SimlistTest185_83.txt

NR=CR190
filename=SimlistTest190_83.txt
#filename=SimlistTest190.txt

#NR=CR190_77deg
#filename=SimlistTest190_77_40_10.txt
#filename=SimlistTest190_77_40_10_small.txt

#NR=CR190_72deg
#filename=SimlistTest190_72_40_10.txt


#NR=CR190_77deg_flat
#filename=SimlistTest190_77_40_0.txt
#filename=SimlistTest190_77_40_0_small.txt

CODE=${NR}_0 #number of your first simulation
CODE2=${NR}_0 #the same number again


FILES=70 #in general 50 protons + 20 iron


################ STARTING THE ANALYSIS ################

  cd ${SCRIPTDIR}

  DATADIR=./CR-Sim/${NR} # the folder of all the simulations
  OTFD=./E${CODE2}_500m/ # Outputfile folder
  


  
#### This has to be done ONE time after the simulations are finished: renames files and put inputfiles in an special folder
 for (( i=0; i < ${FILES}; i++ )); do
      SIMDIR=${DATADIR}/${NR}_${i}

    FILE=${SIMDIR}/a0.trace

    if [ ! -f "$FILE" ]
    then
        echo "File $FILE does not exists"
        python splitZhairesFields_all.py ${SIMDIR}/ ${NR}_${i}
    fi
done

  
  if ! [ -d "$OTFD" ]; then
 	# Control will enter here if $DIRECTORY exists.
 	  mkdir ${OTFD}
 
  fi
  

  MAPDIR=GRAND_antenna.list # List of GRAND antenna positions
  python main.py ${DATADIR} ${CODE2} ${OTFD}/DAT${CODE2}-pickle.dat ${FILES} 50 200 ${OTFD} ${MAPDIR} ${filename}

### particle file not included
### lower freq:50 - high freq.:200  in MHz
