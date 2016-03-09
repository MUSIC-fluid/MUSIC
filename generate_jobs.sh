#!/usr/bin/env bash

workingfolder= /gs/project/cqn-654-ad/mayank/runs
#mkdir $workingfolder

(
  cd /gs/project/cqn-654-ad/mayank/runs

  for ii in {1..5}
  do 
      echo "generating job-$ii..."
	mkdir job-$ii
      cp -r /gs/project/cqn-654-ad/mayank/music job-$ii
      (
          cd job-$ii/music
          ./generate_music_inputfile.py -iev $ii
      )
  done
)
