#!/bin/bash

##########################################################################################
#GENERAL INSTRUCTIONS:                                                                   #
#You should take care of having the following ingredients in order to have this recipe   #
#working: run card and proc card (in a "cards" folder), MadGraph release, this script    #
#all in the same folder!                                                                 #
#Important: Param card is not mandatory for this script                                  #
##########################################################################################

##########################################################################################
#DISCLAIMER:                                                                             #
#This script has been tested in CMSSW_6_2_11 on slc6 and there is no guarantee it should work in   #
#another releases.                                                                       #
#To try in another releases you should adapt it taking into account to put the correct   #
#parameters as the architechture, the release names and compatibility issues as decribed #
#in the comments in this script                                                          #
#Additionally, this script depends on the well behaviour of the lxplus machines, so if   #
#one job fails the full set of jobs will fail and then you have to try again             #
#Any issues should be addressed to: cms-madgraph-support-team_at_cernSPAMNOT.ch          #
##########################################################################################

##########################################################################################
#For runnning, the following command should be used                                      #
#bash create_gridpack_template.sh NAME_OF_PRODCUTION QUEUE_SELECTION                     #
#Or you can make this script as an executable and the launch it with                     #
#chmod +x create_gridpack_template.sh                                                    #
#./create_gridpack_template.sh NAME_OF_PRODCUTION QUEUE_SELECTION                        #
#by NAME_OF_PRODUCTION you should use the names of run and proc card                     #
#for example if the cards are bb_100_250_proc_card_mg5.dat and bb_100_250_run_card.dat   #
#NAME_OF_PRODUCTION should be bb_100_250                                                 #
#for QUEUE_SELECTION is commonly used 1nd, but you can take another choice from bqueues  #
#If QUEUE_SELECTION is omitted, then run on local machine only (using multiple cores)    #
##########################################################################################

#set -o verbose

echo "Starting job on " `date` #Only to display the starting of production date
echo "Running on " `uname -a` #Only to display the machine where the job is running
echo "System release " `cat /etc/redhat-release` #And the system release

#First you need to set couple of settings:

# name of the run
name=${1}

# which queue
queue=$2

#________________________________________
# to be set for user spesific
# Release to be used to define the environment and the compiler needed

#For correct running you should place at least the run and proc card in a folder under the name "cards" in the same folder where you are going to run the script

export PRODHOME=`pwd`
AFSFOLD=${PRODHOME}/${name}
# the folder where the script works, I guess
AFS_GEN_FOLDER=${PRODHOME}/${name}
# where to search for datacards, that have to follow a naming code: 
#   ${name}_proc_card_mg5.dat
#   ${name}_run_card.dat
CARDSDIR=${PRODHOME}/cards
# where to find the madgraph tarred distribution
MGDIR=${PRODHOME}

#MG=MG5_aMC_v2.1.2.tar.gz
MG=MG5_aMC_v2.2.1.tar.gz
MGSOURCE=${PRODHOME}/$MG

SYSCALC=SysCalc_V1.1.0.tar.gz
SYSCALCSOURCE=https://cms-project-generators.web.cern.ch/cms-project-generators/$SYSCALC

MGBASEDIR=MadGraph5_v1_5_14

############################
#Create a workplace to work#
############################
#scram project -n ${name}_gridpack CMSSW ${RELEASE} ; cd ${name}_gridpack ; mkdir -p work ; cd work
WORKDIR=`pwd`
#eval `scram runtime -sh`


#############################################
#Copy, Unzip and Delete the MadGraph tarball#
#############################################
#MGSOURCE=${AFS_GEN_FOLDER}/${MG}

#wget --no-check-certificate ${MGSOURCE}
#cp ${MGSOURCE} .
#wget --no-check-certificate ${MGSOURCE}
#tar xzf $PRODHOME/${MG}
#rm $MG

#############################################
#Apply any necessary patches on top of official release
#############################################

#patch -l -p0 -i $PRODHOME/patches/mgfixes.patch
#patch -l -p0 -i $PRODHOME/patches/models.patch


#mv $MGBASEDIR mgbasedir
#MGBASEDIR=mgbasedir

cd $MGBASEDIR
#cp -pr $PRODHOME/Vector_Triplet_UFO models/.


echo "set auto_update 0" > mgconfigscript
echo "set automatic_html_opening False" >> mgconfigscript

if [ -n "$queue" ]; then
    echo "set run_mode  1" >> mgconfigscript
    echo "set cluster_type lsf" >> mgconfigscript
    echo "set cluster_queue $queue" >> mgconfigscript
    echo "set cluster_status_update 60 30" >> mgconfigscript
    echo "set cluster_nb_retry 3" >> mgconfigscript
    echo "set cluster_retry_wait 300" >> mgconfigscript 
else
    echo "set run_mode 2" >> mgconfigscript
fi

echo "save options" >> mgconfigscript

./bin/mg5 mgconfigscript

#get syscalc and compile
#cp ${SYSCALCSOURCE} .
#wget --no-check-certificate ${SYSCALCSOURCE}
#tar xzf ${SYSCALC}
#rm $SYSCALC

#cd SysCalc
#sed -i "s#INCLUDES =  -I../include#INCLUDES =  -I../include -I${LHAPDFINCLUDES}#g" src/Makefile
#sed -i "s#LIBS = -lLHAPDF#LIBS = ${LHAPDFLIBS}/libLHAPDF.a -lgfortran#g" src/Makefile
#make

#cd $WORKDIR

if [ "$name" == "interactive" ]; then
  exit 0
fi

echo `pwd`


########################
#Locating the proc card#
########################
if [ ! -e $CARDSDIR/${name}_proc_card.dat ]; then
	echo $CARDSDIR/${name}_proc_card.dat " does not exist!"
	#exit 1;
else
	cp $CARDSDIR/${name}_proc_card.dat ${name}_proc_card.dat
fi

########################
#Run the code-generation step to create the process directory
########################

./bin/mg5 ${name}_proc_card.dat



cd $name

#######################
#Locating the run card#
#######################
if [ ! -e $CARDSDIR/${name}_run_card.dat ]; then
echo $CARDSDIR/${name}_run_card.dat " does not exist!"
#exit 1;
else
cp $CARDSDIR/${name}_run_card.dat ./Cards/run_card.dat
fi      

cp $CARDSDIR/${name}_param_card.dat ./Cards/param_card.dat

#######################
#Locating the madspin card#
#######################
domadspin=0
if [ ! -e $CARDSDIR/${name}_madspin_card.dat ]; then
        echo $CARDSDIR/${name}_madspin_card.dat " does not exist! MadSpin will not be run."
else
        cp $CARDSDIR/${name}_madspin_card.dat ./Cards/madspin_card.dat
        domadspin=1
fi


./bin/generate_events -f pilotrun

  
echo "End of job"
  
exit 0

