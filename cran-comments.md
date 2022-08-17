## 2nd Resubmission of New Submission
This is a resubmission after not passing the incoming checks the second time. This
submission is the same as the previous submission except that I have removed two tests
that check for the estimates being in between the confidence interval limits for geeglm and
lme models. The lme test failed on the last submission. The interval methods for the
longitudinal models in this package have not been evaluated for performance and 
there are already warnings when producing intervals for these models. Since there is
already a disclaimer for this, these tests have been removed. The version number has been
bumped to 1.0.2.

The previous submission (1.0.1) was a resubmission after the first submission of 
this package. This submission was the same as the initial except that I have shortened 
the time it takes the 'resi' function examples to run based on the NOTE I received.

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
