# GitTests failed, what now...?!

## Core concepts for beginners, Travis CI

**Continuous Integration** (CI) is the practice of merging in small code changes frequently, rather than merging in a large change at the end of a development cycles. The goal is to build healthier code by developing and testing in smaller increments. *Travis CI* is a platform supporting the CI of projects hosted on GitHub, by automatically building and testing code changes and providing feedback.

**CI builds and automation:** Travis carries out a series of tasks to build and test the code in a brand new virtual environment. The build is considered as passed if none of the tasks fails, and as broken otherwise.

... more at https://docs.travis-ci.com/user/for-beginners/



## What is coveralls.io

Source: [Covaralls documentation](https://docs.coveralls.io).

`coveralls.io` is a web service to help track code coverage over time, and ensure that all new code is fully covered (by unit tests). **The code must be hosted on GitHub** or similar services (BitBucket, GitLab).

**Important remark:** Currently no relevant line of my code is touched, according to `coveralls.io`. This could be due to the use of multi processing. So, for the next push I'll set `use_multiprocessing` to `False` in the configuration parameter file.