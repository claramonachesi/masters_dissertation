#!/bin/bash
# Adaptado de https://github.com/monteirotorres/docking-scripts/blob/main/bash/vbc
echo "#######################################################################################"
echo "#                                  VinaBatchCluster                                   #"
echo "#                                   TORRES, P.H.M.                                    #"
echo "#                                     11-9-2017                                       #"
echo "#######################################################################################"
echo "# This script is intended to perform several docking runs using VINA, then rescore    #"
echo "# the poses using oddt-cli and cluster them using GROMACS gmx cluster. The output is  #"
echo "# a sdf file containing a column corresponding to the cluster ID. Please make sure    #"
echo "#                     oddt, gromacs, babel and vina are installed.                    #"
echo "#######################################################################################"
echo "# Usage:                                                                              #"
echo "# $ ./vcb lig rec vinaconf #runs cutoff                                               #"
echo "# To remove previously generated files:                                               #"
echo "# $ ./vbc rm                                                                          #"
echo "#######################################################################################"
echo


# GET PARAM FROM CMD LINE OR ASK AND CHECK if EXISTS

ligand=${1}
receptor=${2}
config=${3}
runs=${4}
cutoff=${5}

if [ -z ${1+x} ]; then
    echo "Please provide the name of the ligand pdbqt file."
    read ligand
fi

if [ ${1} == 'rm' ]; then
    rm -rf vbcrun* allposes.pdb allposes.pdbqt rmsd-dist.xvg rmsd-clust.xpm clusterinfo.out cluster.log clusters.pdb *clid.xvg RFScore_v3_pdbbind2015.pickle *Rescored.sdf
    echo "Files successfully removed."
    exit
fi

if [ ! -f ${ligand} ]; then
    echo "File ${ligand} not found!"
    exit
fi

if [ -z ${2+x} ]; then
    echo "Please provide the name of the receptor pdbqt file."
    read receptor
fi

if [ ! -f ${receptor} ]; then
    echo "File ${receptor} not found!"
    exit
fi

if [ -z ${3+x} ]; then
    echo "Please provide the name of the vina configuration file."
    read config
fi

if [ ! -f ${config} ]; then
    echo "File ${config} not found!"
    exit
fi

if [ -z ${4+x} ]; then
    echo "How many times should VINA be run with the conf.txt parameters?"
    read runs
fi

if [ -z ${5+x} ]; then
    echo "Choose RMSD cut-off for clustering final poses"
    read cutoff
fi

ligand=`echo $ligand | sed 's/.pdbqt//'`
receptor=`echo $receptor | sed 's/.pdbqt//'`

# DOCK WITH VINA CONVERT OUTPUT USING BABEL
# MAKE MULTI-PDB AND MULTI-PDBQT

for i in `seq 1 ${runs}`; do
    mkdir vbcrun$i
    cd vbcrun$i
    cp ../*.pdbqt .
    cp ../dock_conf.txt .
    echo "Run: ${i}"
    echo ""
    vina --config ${config} 
    obabel -ipdbqt ${ligand}_out.pdbqt -opdb -O ${ligand}_out.pdb
    sed -i "/VINA/a NRUNS ${i}" ${ligand}_out.pdbqt
    sed -i "/VINA/a NRUNS ${i}" ${ligand}_out.pdb
    cat ${ligand}_out.pdb >> ../allposes.pdb
    cat ${ligand}_out.pdbqt >> ../allposes.pdbqt
    cd ..
done

# CLUSTER W/ GMX

gmx cluster -f allposes.pdb -s allposes.pdb -clid ${ligand}-clid.xvg -nofit -method gromos -cutoff ${cutoff} -cl clusters.pdb <<EOF
1
1
EOF
tail -n +21 ${ligand}-clid.xvg > clusterinfo.out

# ADD CLUSTER AS REMARK IN PDB

sed -i '/VINA/a inserirclusteraqui' allposes.pdbqt
cat clusterinfo.out | while read a b; do 
    sed -i "0,/inserirclusteraqui/s//CLUSTR\ $b/" allposes.pdbqt  ;
done
 # RESCORE W/ ODDT

oddt_cli -i pdbqt --score rfscore_v3_pdbbind2015 --receptor ${receptor}.pdbqt -o sdf -O ${ligand}-Rescored.sdf allposes.pdbqt

