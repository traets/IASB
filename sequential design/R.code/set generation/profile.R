##profile

Rprof("profile1.out", line.profiling=TRUE)
eval(parse(file = "C:/Users/u0105757/Desktop/sequential_design.sas/R.code/IASB_R/DBerror_MODfederov.9.0.R", keep.source=TRUE))
Rprof(NULL)

summaryRprof("profile1.out", lines = "show")
