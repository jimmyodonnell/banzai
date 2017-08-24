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
		install.packages(package, repos = 'https://cloud.r-project.org/')
	}
}

lib <- .libPaths()

if(any(! package_req %in% installed.packages())){
	pkg_missing <- package_req[(! package_req %in% installed.packages())]
	long_msg <- paste("Please check your internet connection and try again.", 
	"If you are certain you have installed the required packages, they may be ", 
	"in a different location. You could install them by hand and try again; ", 
	"just be sure you install them into this directory: ")
	msg <- paste(c(
		"Unable to find or install required R packages:", "",  
	  pkg_missing, "", strwrap(long_msg), "", lib[1], "", 
		"You can accomplish this in R using:", "", 
		"install.packages('ThePackageName', lib = 'paste_the_above_path_here')", 
		"\n"
	), collapse = "\n")
	stop(msg)
}
