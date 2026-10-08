# SCODE

R implementation of SCODE for regulatory network inference from time-course expression data.

The benchmark adapter in `code/methods/SCODE.R` calls `run_scode()` from `SCODE.R`
in a separate R process. Settings are defined in the adapter; the vendored implementation
uses MASS and exports `A.txt`, `B.txt`, `W.txt` and `RSS.txt`.

[Upstream source](https://github.com/hmatsu1226/SCODE) ·
[Method paper](https://academic.oup.com/bioinformatics/article/33/15/2314/3100331)

Copyright (c) 2016 Hirotaka Matsumoto. Released under the MIT license.
