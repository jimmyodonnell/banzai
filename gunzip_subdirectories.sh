#!/bin/bash

my_dir="$(dirname "$0")"
find "${my_dir%/*}" -type f -name '*.gz' -exec gunzip "{}" \;
