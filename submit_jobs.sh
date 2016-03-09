#!/usr/bin/env bash


(
 cd /gs/project/cqn-654-ad/mayank/runs/
      
   for ii in `ls | grep "job-"`
   do
       echo "submitting $ii ..."
       (
           cd $ii/music
           qsub submit_full_job.pbs
       )
   done
)
