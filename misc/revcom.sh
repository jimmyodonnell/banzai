#!/usr/bin/env bash
# a simple function to reverse complement DNA strings
# Jimmy O'Donnell <jodonnellbio@gmail.com>

revcom (){

  echo $1 |\
    rev   |\
    tr ACGTWSMKRYBDHVNacgtwsmkrybdhvn TGCASWKMYRVHDBNtgcaswkmyrvhdbn

}

# alt:
#  tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]"
