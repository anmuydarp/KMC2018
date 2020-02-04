#!/bin/bash

END=40
thickness=6nm
temp=300K
path=/home/pmuralid/KMC2018/ForJournal/Thickness/${temp}/${thickness}/Code

for i in $(seq 1 $END)
do
    cd $i.data/

    File=${thickness}_${i}_data.txt

    if [ -f $File ]
    then        
        #PF
        A=$( echo "1p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c11- )
        echo -e $A1 >> ${path}/${thickness}_PF.txt 

        #D2E        
        A=$( echo "2p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c7- )
        echo -e $A1 >> ${path}/${thickness}_D2E.txt 

        #Thermionic
        A=$( echo "3p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c23- )
        echo -e $A1 >> ${path}/${thickness}_Therm.txt 
        
        #0 phonon
        A=$( echo "4p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_0phonon.txt 
        
        #1 phonon
        A=$( echo "5p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_1phonon.txt 
        
        #2 phonon
        A=$( echo "6p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_2phonon.txt 

        #3 phonon
        A=$( echo "7p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_3phonon.txt 

        #4 phonon
        A=$( echo "8p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_4phonon.txt 

        #5 phonon
        A=$( echo "9p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_5phonon.txt 

        #Avg Time
        A=$( echo "11p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c14- )
        echo -e $A1 >> ${path}/${thickness}_energy.txt     
        
        #Avg tau
        A=$( echo "12p" | ed -s "$File" )        
        A1=$( echo "$A" | cut -c12- )
        echo -e $A1 >> ${path}/${thickness}_tau.txt 
    else
        echo "File NOT FOUND"
    fi
    cd ${path}
done
echo "Mined"
