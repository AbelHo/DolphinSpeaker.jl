# End-to-end benchmark on simulated pool recordings

```
55 scenarios x 4 s, A0 9000, threshold 490.4; 22440 in-train detections, 99.5% matched to an arrival, 61.1% of them direct-path
simulated vs real: followers 42% (real 51%); median ppeak/threshold 2.11 (real 1.79); SIC >= 2 candidates 15% (real 78%)

selected direction within 3 deg of the emitting click's direct path (%), and the real-data proxies on the simulation
win   tracker       rule                                     all  direct    echo med err |   burst  2-10ms 10-30ms  frame
sic   -             oracle (best candidate)                 65.6    99.9    18.1     0.4 |     NaN     NaN     NaN    NaN
sic   birth 0.1     track lag + followers inherit           80.4    92.0    64.2     0.3 |     0.4     0.5     0.5   64.6
sic   birth 0.1     track lag                               65.0    99.7    16.9     0.4 |     0.4     0.5     0.5   63.1
sic   birth 0.3     track lag + followers inherit           80.3    92.0    64.1     0.3 |     0.4     0.5     0.5   64.6
sic   birth 0.3     track lag                               65.0    99.7    16.8     0.4 |     0.4     0.5     0.5   63.1
sic   default       track lag + followers inherit           80.4    92.0    64.2     0.3 |     0.4     0.5     0.5   64.6
sic   default       longest track + followers inherit       80.2    92.1    63.6     0.3 |     0.4     0.5     0.5   64.9
sic   default       earliest + followers inherit            80.0    91.8    63.8     0.3 |     0.4     0.5     0.5   64.8
sic   default       rank 1 + followers inherit              79.2    91.5    62.2     0.3 |     0.4     0.5     0.5   64.9
sic   default       longest track                           65.2    99.8    17.2     0.4 |     0.4     0.5     0.5   63.4
sic   default       track lag                               65.0    99.7    16.9     0.4 |     0.4     0.5     0.5   63.1
sic   default       earliest                                64.7    99.4    16.5     0.4 |     0.4     0.5     0.5   63.2
sic   default       rank 1                                  64.4    99.7    15.6     0.4 |     0.4     0.5     0.5   63.2
sic   delete 0.5 s  track lag + followers inherit           80.4    92.0    64.2     0.3 |     0.4     0.5     0.5   64.6
sic   delete 0.5 s  track lag                               65.0    99.7    16.9     0.4 |     0.4     0.5     0.5   63.1
sic   gate 4        track lag + followers inherit           80.3    92.0    64.1     0.3 |     0.4     0.5     0.5   64.6
sic   gate 4        track lag                               65.0    99.7    16.8     0.4 |     0.4     0.5     0.5   63.2
sic   gate 9        track lag + followers inherit           80.4    92.0    64.2     0.3 |     0.4     0.5     0.5   64.7
sic   gate 9        track lag                               65.0    99.7    16.9     0.4 |     0.4     0.5     0.5   63.1
sic   q 10          track lag + followers inherit           80.3    92.0    64.2     0.3 |     0.4     0.5     0.5   64.6
sic   q 10          track lag                               65.0    99.7    16.8     0.4 |     0.4     0.5     0.5   63.1
sic   q 40          track lag + followers inherit           80.4    92.0    64.2     0.3 |     0.4     0.5     0.5   64.5
sic   q 40          track lag                               65.0    99.7    16.9     0.4 |     0.4     0.5     0.5   63.1
sicw  -             oracle (best candidate)                 71.4    99.9    31.8     0.4 |     NaN     NaN     NaN    NaN
sicw  birth 0.1     track lag + followers inherit           80.2    91.2    64.9     0.3 |     0.4     0.4     0.5   64.7
sicw  birth 0.1     track lag                               68.0    98.0    26.3     0.4 |     0.4     0.4     0.5   65.3
sicw  birth 0.3     track lag + followers inherit           80.2    91.2    64.8     0.3 |     0.4     0.4     0.5   64.6
sicw  birth 0.3     track lag                               67.9    98.0    26.2     0.4 |     0.4     0.4     0.5   65.2
sicw  default       earliest + followers inherit            80.5    91.1    65.7     0.3 |     0.4     0.4     0.5   64.6
sicw  default       longest track + followers inherit       80.4    91.6    65.0     0.3 |     0.4     0.4     0.5   64.5
sicw  default       track lag + followers inherit           80.2    91.2    64.9     0.3 |     0.4     0.4     0.5   64.7
sicw  default       rank 1 + followers inherit              79.6    91.5    63.0     0.3 |     0.4     0.4     0.5   64.9
sicw  default       longest track                           70.3    98.5    31.0     0.4 |     0.4     0.4     0.5   64.4
sicw  default       rank 1                                  68.6    99.6    25.7     0.4 |     0.4     0.4     0.5   65.0
sicw  default       track lag                               67.9    98.0    26.3     0.4 |     0.4     0.4     0.5   65.3
sicw  default       earliest                                67.4    97.8    25.2     0.4 |     0.4     0.4     0.5   65.6
sicw  delete 0.5 s  track lag + followers inherit           80.2    91.2    64.9     0.3 |     0.4     0.4     0.5   64.7
sicw  delete 0.5 s  track lag                               67.9    98.0    26.3     0.4 |     0.4     0.4     0.5   65.3
sicw  gate 4        track lag + followers inherit           80.2    91.2    64.9     0.3 |     0.4     0.4     0.5   64.6
sicw  gate 4        track lag                               68.0    98.0    26.4     0.4 |     0.4     0.4     0.5   65.2
sicw  gate 9        track lag + followers inherit           80.2    91.2    64.8     0.3 |     0.4     0.4     0.5   64.7
sicw  gate 9        track lag                               68.0    98.0    26.5     0.4 |     0.4     0.4     0.5   65.3
sicw  q 10          track lag + followers inherit           80.2    91.2    64.9     0.3 |     0.4     0.4     0.5   64.7
sicw  q 10          track lag                               68.0    98.0    26.3     0.4 |     0.4     0.4     0.5   65.3
sicw  q 40          track lag + followers inherit           80.2    91.2    64.8     0.3 |     0.4     0.4     0.5   64.6
sicw  q 40          track lag                               67.9    98.0    26.2     0.4 |     0.4     0.4     0.5   65.2

sic track rule echo label: precision 66.3%, recall 25.2% (image detections: 41.9%)
sicw track rule echo label: precision 43.8%, recall 24.2% (image detections: 41.9%)

rank correlation over the rules (default tracker) between accuracy and each proxy (proxy: lower change / higher frame share = better):
  sic   burst +0.14, 2-10 ms -0.36, 10-30 ms -0.29, in frame +0.76 (n rules 8)
  sicw  burst -0.05, 2-10 ms -0.05, 10-30 ms -0.24, in frame -0.74 (n rules 8)
```
