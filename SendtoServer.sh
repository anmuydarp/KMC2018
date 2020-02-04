#!/bin/bash

File="${runfile}.sh"

cat<<EOT>>$File
#!/bin/bash

#SBATCH -n 1 -t 0-4:00
#SBATCH --qos=normal

module load gcc/7.3.0
                                                                                                                                                                       
${path}

EOT
