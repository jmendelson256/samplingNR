# samplingNR 0.2.9002 - dev

* In `opt_nh_nonresp`, moved ... to be before optional args
* In process of adding handling of different types of objectives. Stuff added so far:
- Added description of intended handling, but haven't checked rest of description to make sure advice is good on which combos of args to use.
- Programmed extensive error handling; moved to separate internal fns to try to keep somewhat abstract.  Added extensive testthat functionality; confident it works.
- Programmed new method -- fixed precision. Seems to work but could add some unit tests.
* Tried test cases for new method.  Works so far:
- Took vignette example, calculated variance, then set that as precision target (both ways); same result. Looked over in Excel somewhat.
- Under complete response, can replicate PracTools vignette example.
- Under nonresponse, changes allocation in expected way.
* Outstanding items.
- Re-read documentation -- still correct?
- Add unit testing for fixed precision method.
- Add fn to calculate variance for finite pop mean. Take a look at jmpaper3::calc_EVar_under_nonresp() but may want to remove some of the options that won't be commonly used to avoid user confusion.

# samplingNR 0.2.1

* Numerous documentation edits

# samplingNR 0.2.0

* Ported core functionality from development version (jmpaper3_0.1.0.9016.tar.gz)
* Renamed functions, cleaned up documentation and examples, and added unit tests
* Added in `pevs_adm_2016_rrs` data set
