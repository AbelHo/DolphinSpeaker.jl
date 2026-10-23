# Round-4 benchmark: realistic simulation, inheritance, n dolphins

```
realistic simulation: A0 6000, hydrophone scatter 0.70, diffuse reverb 0.10, off-axis low-pass on; window -110:50; SIC min_rel 0.7

=== 2 dolphins: 27 scenarios x 4 s, 9872 detections; followers 37% (real 51%), >= 2 SIC candidates 70% (real 79%), closure <= 1 45% (real 35%)
    matched detections: direct path 71%, specular echo 29%, diffuse 0%
  direction within 3 deg of the emitting click's direct path (%):
  method                                   all   direct    echo other-dol.  median°
  longest track + inherit all             36.0     41.7    23.9       28.4     18.5
  longest track + smart inherit           35.2     42.0    20.9       40.2     15.1
  track lag                               32.0     44.0     6.6       45.9     18.8
  longest track                           31.1     42.8     6.4       45.4     18.5
  earliest                                30.5     42.2     6.0       43.8     19.6
  rank 1 + smart inherit                  27.0     33.9    12.5       35.6     24.1
  rank 1 + inherit all                    26.7     32.9    13.6       25.3     25.5
  rank 1                                  24.5     34.5     3.3       38.1     25.8
  smart inherit roles: echo label precision 71.9% / recall 35.0% vs specular+diffuse echo detections; 'own' on other-dolphin followers 9.8%
  source separation (pairs of detections <= 0.2 s apart):
  method                       precision  recall     F1    ARI   n true  n found unassigned%
  tracker, n given                 0.966   0.691  0.777  0.331      1.7      1.7        16.3
  tracker, n estimated             0.968   0.701  0.787  0.352      1.7      5.2         4.9
  k-means baseline (k = n)         0.934   0.804  0.851  0.187      1.7      2.0         0.0
  rules: accuracy vs the real-data proxies computed on the simulation (median over scenarios):
  longest track    acc  39.2% | burst   2.2°  2-10 ms   1.6°  10-30 ms   3.6°
  track lag        acc  38.8% | burst   3.6°  2-10 ms   1.6°  10-30 ms   3.6°
  earliest         acc  37.4% | burst  10.8°  2-10 ms  12.3°  10-30 ms  14.8°
  rank 1           acc  34.2% | burst  11.9°  2-10 ms   2.6°  10-30 ms   3.6°

=== 3 dolphins: 14 scenarios x 4 s, 9244 detections; followers 36% (real 51%), >= 2 SIC candidates 81% (real 79%), closure <= 1 39% (real 35%)
    matched detections: direct path 74%, specular echo 26%, diffuse 0%
  direction within 3 deg of the emitting click's direct path (%):
  method                                   all   direct    echo other-dol.  median°
  longest track + smart inherit           53.6     62.2    29.3       39.3      1.8
  longest track + inherit all             53.4     60.4    33.7       23.2      2.0
  track lag                               50.9     66.1     8.2       55.5      2.6
  earliest                                49.1     65.0     4.3       51.7      3.6
  rank 1 + inherit all                    48.9     55.6    30.2       24.8      3.6
  rank 1 + smart inherit                  48.6     56.9    25.2       36.1      3.6
  longest track                           48.5     63.1     7.2       50.9      3.5
  rank 1                                  44.2     58.3     4.5       52.0      9.6
  smart inherit roles: echo label precision 67.3% / recall 35.7% vs specular+diffuse echo detections; 'own' on other-dolphin followers 20.7%
  source separation (pairs of detections <= 0.2 s apart):
  method                       precision  recall     F1    ARI   n true  n found unassigned%
  tracker, n given                 0.928   0.697  0.769  0.355      2.8      2.8        17.6
  tracker, n estimated             0.930   0.731  0.797  0.491      2.8      7.8         3.3
  k-means baseline (k = n)         0.893   0.754  0.801  0.306      2.8      3.0         0.0
  rules: accuracy vs the real-data proxies computed on the simulation (median over scenarios):
  track lag        acc  52.4% | burst   2.3°  2-10 ms   4.7°  10-30 ms   2.0°
  earliest         acc  51.2% | burst   4.2°  2-10 ms  21.4°  10-30 ms  17.6°
  longest track    acc  50.3% | burst   1.5°  2-10 ms   3.0°  10-30 ms   3.7°
  rank 1           acc  44.3% | burst   1.2°  2-10 ms  14.1°  10-30 ms  15.5°
```
