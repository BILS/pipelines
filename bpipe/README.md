#		Bpipe pipeline collection
#		########################


 FOR NBIS's employees full details are available under Infrastructure>Pipelines>Bpipe on the Redmine
-----------------------------------------------------------------------------------------------------


[1. For the impatient](#for-the-impatient)</br>
[2. Disclaimer](#disclaimer)</br>
[3. Installation](#installation)</br>
[4. Usage information](#usage-information)</br>
[5. Licence](#licence)</br>
[6. Acknowledgment](#acknowledgment)

</br>
##1. For the impatient

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

##2. Disclaimer

In the case where this work is used for communications (presentation, publication, etc.), we invite you to acknowledge BILS within your communication according to this example: "Support by BILS (Bioinformatics Infrastructure for Life Sciences) is gratefully acknowledged."

##3. Installation

Bpipe does not need to be installed.  However, in order to run the
various pipelines, you will need to make sure that all the dependencies
are available and configured in the respective pipeline.config files
(more below).  In addition, you will need to make sure that the location
of the pipeline modules is exported to $BPIPE_LIB.  For a quick way to
find out how, run the script 'setup.sh' that is included with this code.

Also make sure that you generate the config files for your pipeline. 

bpipe.config.template

pipeline.config.template

To generate config files, you can run:

$PIPELINE_GIT/config/bpipe_config -p Name_of_Pipeline -c

(to find out which pipelines are available, run bpipe_config -l).

$BPIPE_GIT/config/bpipe_config -b

###3.1 Bpipe.config

This file enables the translation of a pipeline to a resource manager
or job queue.  Common examples include Slurm, LSF or SGE.  The template
includes some suggestions on parameters for Slurm and LSF.

You can then define the required resources for each program that you are
going to run in the pipeline.  For more details, please refer to the
bpipe documentation on www.bpipe.org.

###3.2 pipeline.config

The pipeline config file holds variables that specify the location of
binaries, input files or values required by the various pipeline stages.
Each variable should be documented and it is for you to fill it out
correctly.

##4. Usage information

In order to start a pipeline run, type:

    bpipe run path/to/pipeline.bpipe

If you do not specify one or several input files, the pipeline will
inform you about the expected input data.  With that information, do:

    bpipe run path/to/pipeline.bpipe input

For an RNA-seq workflow, this could be:

    bpipe run $GIT/pipelines/pipelines/RNAseqAlignQuantify.bpipe *.fq.gz

Where *.fq.gz specifies a pattern of gzipped fastq files.  The pipeline
will separate these into individual branches and run them in parallel
where applicable.

##5. Licence
This project is licensed under the terms of the GNU GPLv3 license. Please read the [LICENSE](LICENSE) file for more information.

##6. Acknowledgment
We acknowledge the bpipe developers => [Bpipe project](https://github.com/ssadedin/bpipe).</br>
Sadedin, S. P., Pope, B. & Oshlack, A. Bpipe: A tool for running and managing bioinformatics pipelines. Bioinformatics 28, 1525–1526 (2012).</br>
Special thanks to Marc Hoeppner whose start this bpipe pipeline library.
