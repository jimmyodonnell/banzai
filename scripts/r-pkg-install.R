#!/usr/bin/env Rscript

arguments <- commandArgs(TRUE)

pkg_file <- arguments[1]

################################################################################
# LOAD PACKAGES
################################################################################
package_req <- read.table(pkg_file, stringsAsFactors = FALSE)[,1]
# res <- vector()
for(package in package_req){
	if(! package %in% installed.packages()){
		install.packages(package)
	}
}
if(any(! package_req %in% installed.packages())){
	pkg_missing <- package_req[(! package_req %in% installed.packages())]	
	msg <- paste0(c(
		"Unable to find or install required R packages:", "\n", 
	  "\t", pkg_missing, "\n", 
	  "Please check your internet connection and try again."
	))
	stop(msg)
}
