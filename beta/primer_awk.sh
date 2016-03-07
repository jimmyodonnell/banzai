#!/usr/bin/env bash


 # grep for AAA and BBB and CCC (in any order on the same line)
 awk '/AAA/ && /BBB/ && /CCC/'

 # grep for AAA and BBB and CCC (in that order)
 awk '/AAA.*BBB.*CCC/'
 
 
  # grep for AAA and BBB and CCC (in that order)
 awk '/.{0,range}primer1.*primer1rc.{0,range}/'
 
 