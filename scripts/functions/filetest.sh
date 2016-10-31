#!/usr/bin/env bash

filetest () {
  if [[ -s "${1}" ]]; then
    timestamp "The following file is empty or absent:"
    echo "${1}"
    echo "This is a critical file, so the script will exit."
    exit
  fi
}
