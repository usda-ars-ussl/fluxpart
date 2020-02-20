# v0.2.8 (2020-02-20)
- Fix bug affecting integrity check of TOB datafiles

# v0.2.7 (2019-12-27)
- Avoid divide-by-zero when Fcr/Fcp = -1
- Update for Python 3.7 changes in StopIteration exceptions raised by generators

# v0.2.6 (2019-11-27)
- Add "opt" model for intercellular leaf CO2
- Additional checks for missing high frequency data
- Additional outputs added to high freqency data summary
- Improved error messaging

# v0.2.5 (2019-10-31)
- Add support for q and c data being molar ratio instead of mass concentration
- Fix typo in docs figure legend (GPP -> NEE)

# v0.2.4 (2019-04-30)
- No longer saving/pickling fluxpart meta data, was not working with pandas 0.24  

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
