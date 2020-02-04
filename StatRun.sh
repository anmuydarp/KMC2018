#!/bin/bash

END=40
thickness=6nm
temp=300K

for i in $(seq 1 $END)
do
    mkdir $i.data

    ipfile=KMCmain_$i
    path=/home/pmuralid/KMC2018/ForJournal/Thickness/${temp}/${thickness}/Code/$i.data/KineticMC
    runfile=${thickness}_$i

    echo $ipfile
    export ipfile
    export i
    export path
    export runfile

    sh ./CreateKMCMain.sh   
    sh ./SendtoServer.sh
    
    mv ${ipfile}.cpp ${runfile}.sh $i.data/
    cp -r Math.cpp Initialize.cpp CreateDefects.cpp Device.cpp TransitionTable.cpp TransitionRateCalc.cpp TransitionRateFunctions.cpp $i.data/
    cp KMCLoop.cpp KMCEvents.cpp Write_Calc.cpp Lifetime.cpp $i.data/
    cp *.h Makefile $i.data/

    cd $i.data/

    mv ${ipfile}.cpp KMCmain.cpp

    make

    pwd

    sbatch ${runfile}.sh

    cd /home/pmuralid/KMC2018/ForJournal/Thickness/${temp}/${thickness}/Code

done

