# metacom_dis_fit
Code for Dispersal mode and spatial heterogeneity shape the interaction between adaptation and dispersal in multitrophic metacommunities

This is a very rough version, check GIT for a more polished one, but this is what we run during research 

We run the code on a cluster based on
OS: Debian GNU/Linux 12 (bookworm) x86_64
Mostly using Xeon(R) Gold 6140 CPU and Xeon(R) Silver 4210R CPUS (about 250 cores, with at least 4GB per core )
Queuing system: HTcondor 25.4.0 

We tested the code with Matlab R2016a, but it should be compatible with many other versions, as we did not use unusual features.

Code usage:
The code leaves logs in a ./log directory and the results in a ./results_paper directory.

To run the code use

condor_submit <jobfile>

The results use the job cluster name as part of the number, so you have to write that number down.


The jobfiles for the main article are

jobs_ns25_ne50_Tr0.5_beta0.99_SI.job
jobs_ns25_ne50_Tr0.5_beta0_SI.job
jobs_ns25_ne50_Tr15_beta0.99_SI.job
jobs_ns25_ne50_Tr15_beta0_SI.job
jobs_ns25_ne50_Tr5_beta0.99_SI.job
jobs_ns25_ne50_Tr5_beta0_SI.job

(in spite of the SI suffix)

After the jobs finish, the plots are made using a post processing program (no need to use a cluster). See the post_proc directory
