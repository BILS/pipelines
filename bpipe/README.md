#		Bpipe pipeline collection
#		########################


 FOR NBIS's employees full details are available under Infrastructure>Pipelines>Bpipe on the Redmine
-----------------------------------------------------------------------------------------------------


[1. For the impatient](#1-for-the-impatient)</br>
[2. Disclaimer](#2-disclaimer)</br>
[3. Installation](#3-installation)</br>
&nbsp;&nbsp;&nbsp;[3a Bpipe config file](#3a-bpipe-config-file)</br>
&nbsp;&nbsp;&nbsp;[3b Pipeline config file](#3b-pipeline-config-file)</br>
[4. Usage information](#4-usage-information)</br>
[5. Licence](#5-licence)</br>
[6. Acknowledgment](#6-acknowledgment)</br>

## 1 For the impatient

This is for those who have a basic idea of how bpipe works:

    -   Make sure that the module directory of this code base is
        exported to $BPIPE_LIB (run ./setup.sh to check).

    -   Generate a pipeline config and bpipe config file using the 
    	included tool (config/bpipe_config)

    -   Rename the files to bpipe.config and pipeline.config.

    -   Fill out the files.

    -   Make sure that the various tools are loaded into $PATH (unless
        you specify full paths to them in the config file).

    -   Run pipeline on your input data.

## 2 Disclaimer

In the case where this work is used for communications (presentation, publication, etc.), we invite you to acknowledge NBIS within your communication according to this example: "Support by NBIS (National Bioinformatics Infrastructure Sweden) is gratefully acknowledged."

## 3 Installation

### General information
Bpipe does not need to be installed.  However, in order to run the
various pipelines, you will need to make sure that all the dependencies
are available and configured in the respective pipeline.config files
(more below).  In addition, you will need to make sure that the location
of the pipeline modules is exported to $BPIPE_LIB.  For a quick way to
find out how, run the script 'setup.sh' that is included with this code.

Also make sure that you generate the config files for your pipeline: __bpipe.config.template__ and __pipeline.config.template__.

### Step by step:
 * **clone the repo and move in it:**
 
       git clone https://github.com/NBISweden/pipelines.git
       cd pipelines/bpipe/
       
 * **check BPIPE variables (to make everything smoother):**
 
       ./setup.sh
       
 * **accoding to the ouput from the ` ./setup.sh` add the missing variables in your bash_profile and source it:**
 
       source ~/.bash_profile
 
 * **Bpipe config file**

This file (bpipe.config) enables the translation of a pipeline to a resource manager or job queue. 
Common examples include Slurm, LSF or SGE.  The template includes some suggestions on parameters for Slurm and LSF.

You can then define the required resources for each program that you are
going to run in the pipeline.  For more details, please refer to the
bpipe documentation on www.bpipe.org.

To generate the bpipe config file run:  

    $BPIPE_CONFIG/bpipe_config -b

 * **Pipeline config file**

The pipeline config (pipeline.config) file holds variables that specify the location of
binaries, input files or values required by the various pipeline stages.
Each variable should be documented and it is for you to fill it out
correctly.

To generate the pipeline config file, run:

    $BPIPE_CONFIG/bpipe_config -p Name_of_Pipeline -c

To look at available pipelines, run :

    $BPIPE_CONFIG/bpipe_config –l


## 4 Usage information

In order to start a pipeline run, type:

    bpipe $BPIPE_PIPELINE/pipeline.bpipe

If you do not specify one or several input files, the pipeline will
inform you about the expected input data.  With that information, do:

    bpipe run $BPIPE_PIPELINE/pipeline.bpipe input

For an RNA-seq workflow, this could be:

    bpipe run $BPIPE_PIPELINE/RNAseqAlignQuantify.bpipe *.fq.gz

Where *.fq.gz specifies a pattern of gzipped fastq files.  The pipeline
will separate these into individual branches and run them in parallel
where applicable.

## 5 Licence
This project is licensed under the terms of the GNU GPLv3 license. Please read the [LICENSE](LICENSE) file for more information.

## 6 Acknowledgment
We acknowledge the bpipe developers => [Bpipe project](https://github.com/ssadedin/bpipe).</br>
Sadedin, S. P., Pope, B. & Oshlack, A. Bpipe: A tool for running and managing bioinformatics pipelines. Bioinformatics 28, 1525–1526 (2012).</br>
Special thanks to Marc Hoeppner whose start this bpipe pipeline library.
