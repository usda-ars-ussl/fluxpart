# v0.2.3 (2019-04-30)
- Fix save/read of fluxpart results for pandas 0.24  

# v0.2.2 (2019-04-25)
- Fix filepath bug occuring when datafiles are GHG format

# v0.2.1 (2019-03-18)
- Add support for date and time being in separate datafile columns
- Add suport for GHG datafile format
- Improved root finding 

# v0.2.0 (2018-09-19)
- New primary module interface is ``fvs_partition``
- Default csv data column ordering changed
- Datafiles can be specified using file/path wildcards, w/ optional time-sorting
- Flux partioning interval indendent of datafile time intervals 
- Specify plant & tower heights or measured wue via file or callable
- Add option to make all fluxes nonstomatal after/before sunset/sunrise
- Add support for TOB datafile format
- New format for results, stored in pandas dataframe, using common units

# v0.1.0 (2018-05-05)
- Initial release
