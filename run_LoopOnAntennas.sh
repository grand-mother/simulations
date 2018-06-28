#!/bin/sh
#. /home/renault/.bash_profile
#module load intel/17.0-python-2.7.12

script_folder="/Users/nrenault/Desktop/GRAND/scripts_clean/"
date

#############################################################################################
##                                                                                         ##
##        Script written to run zhaires simulations, antenna response computations,        ##
##                   filter signals and perform some 2D map plotting                       ##
##                                                                                         ##                                                                     \
#############################################################################################

# For any zhaires simulations (CR, ToyModel, ...), only need to provide the path to the directory containing the simulation directory of each shower.
if [ $1 = "else" ]; then
  #Set directories
  home_folder="/Users/nrenault/Desktop/GRAND/"
  output_folder=$2 #Need to input the full path, not just a part of it starting in current directory
  echo "output_folder = "$output_folder
  cd $output_folder

  # Run raw simulation data file splitting, voltage computing, signal filtering, 2D map plotting, detection check
  echo "Computing voltage for each antenna"
  lfreq = 50 #MHz
  hfreq = 200 #MHz
  shower_files=${output_folder}"/"*"/"
  for shower_input in $shower_files; do
    shower_timefresnel=${shower_input}"/timefresnel-root.dat"
    if [ -f $shower_timefresnel ]; then

      #Comment/Uncomment the following line to run or not the timefresnel-root.dat file splitting
      python ${script_folder}/splitZhairesFields_all.py ${shower_input}
      echo "Compute voltage..."
      inp_file=${shower_input}/*.inp

      #Comment/Uncomment the following lines to run or not the voltage computing  and Efield, and voltage 2D map plotting.
      python ${script_folder}/computevoltage_HorAnt.py txt 1 ${shower_input}/split/ ${inp_file}
      python ${script_folder}/Efield_2Dmap.py ${shower_input}"/" 0 0
      python ${script_folder}/voltage_2Dmap.py ${shower_input}"/" 0 0

      #Comment/Uncomment the following lines to run or not the frequency filtering of the voltage (or the Efield with the 'Efield' keyword)
      trace_file=${shower_input}/split/*.trace
      for antenna_input in ${trace_file}; do
        python ${script_folder}/processSim_Nicolas.py ${shower_input}/split/ ${antenna_input} $lfreq $hfreq 'voltage'
      done

      #Comment/Uncomment the following line to run or not the plotting of the filtered voltage
      python ${script_folder}/voltage_2Dmap.py ${shower_input}"/" $lfreq $hfreq
    fi
  done

  #Comment/Uncomment the following lines to perform the detection results for each shower
  python ${script_folder}/check_detection.py toymodel #CR and GP300 is an other keyword to use for CR 300km2 flat or CR GP300 simulations
  #python ${script_folder}/check_detection.py toymodel $lfreq $hfreq

############################################################################
############################################################################
############################################################################

# Dedicated to RadioMorphing simulation, only need to provide the path to the directory containing the simulation directory of each shower, and the path to the directory containing the corresponding zhaires input files.
#This script part was dedicated to the for RM and ZHAireS simulations  within the idea to compare them
elif [ $1 = "RM" ]; then
  #Set directories
  home_folder="/Users/nrenault/Desktop/GRAND/"
  work_folder=$home_folder"RadioMorphing/simus_RM/"
  output_folder=$2
  inp_folder=$3
  echo "output_folder = "$output_folder
  cd $output_folder

  # Run raw simulation data file splitting, voltage computing, signal filtering, 2D map plotting, detection check
  echo "Computing voltage for each antenna"
  lfreq = 50 #MHz
  hfreq = 200 #MHz
  shower_files=${output_folder}"/"*"/"
  for shower_input in $shower_files; do
    shower_number=$(basename ${shower_input} .inp)
    echo "Compute voltage for shower number: "${shower_number}

    #Comment/Uncomment the following lines to run or not the voltage computing  and Efield, and voltage 2D map plotting.
    inp_file=$inp_folder/*-${shower_number}.inp #$inp_folder"/D40000m-Z10deg-"${shower_number}".inp"
    python ${script_folder}/computevoltage_HorAnt.py txt 1 ${shower_input} ${inp_file}
    python ${script_folder}/Efield_2Dmap.py ${shower_input}"/" 0 0
    python ${script_folder}/voltage_2Dmap.py ${shower_input}"/" 0 0

    #Comment/Uncomment the following lines to run or not the frequency filtering of the voltage (or the Efield with the 'Efield' keyword)
    trace_file=${shower_input}/split/*.trace
    for antenna_input in ${trace_file}; do
      python ${script_folder}/processSim_Nicolas.py ${shower_input}/split/ ${antenna_input} $lfreq $hfreq 'voltage'
    done

    #Comment/Uncomment the following line to run or not the plotting of the filtered voltage
    python ${script_folder}/voltage_2Dmap.py ${shower_input}"/" $lfreq $hfreq
  done

  #Comment/Uncomment the following lines to perform the detection results for each shower
  python ${script_folder}/check_detection.py RM
  #python ${script_folder}/check_detection.py RM" $lfreq $hfreq
fi

############################################################################
echo "Pipeline finished"
cd ..
date
