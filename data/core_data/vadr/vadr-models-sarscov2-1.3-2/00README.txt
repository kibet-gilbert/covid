November 2021
vadr-models-corona-1.3-2

https://github.com/ncbi/vadr

Instructions for SARS-CoV-2 annotation using VADR:
https://github.com/ncbi/vadr/wiki/Coronavirus-annotation

VADR documentation:
https://github.com/ncbi/vadr/blob/master/README.md

See RELEASE-NOTES.txt for details on changes between model versions. 

------------

Recommended command for SARS-CoV-2 annotation 
(used by GenBank as of writing with vadr 1.3):

v-annotate.pl \
--split --cpu 8 \
--glsearch \
-s -r --nomisc \ 
--mkey sarscov2 \
--lowsim5seq 6 --lowsim3seq 6 \
--alt_fail lowscore,insertnn,deletinn \
--mdir <PATH-TO-THIS-MODEL-DIR> \
<fastafile> <outputdir>

The '--split --cpu 8' options will run v-annotate.pl multi-threaded on
8 threads. To change to '<n>' threads use '--split --cpu <n>', but
make sure you have <n> * 2G total RAM available. 

To run single threaded remove the '--split --cpu 8' options.

------------

Contact eric.nawrocki@nih.gov for help.
