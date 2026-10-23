# Detection timing diagnostics

Source: `/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device`. Directions: Impulse6 (SRP-PHAT) and Impulse5 (GCC-PHAT guided) pixels inverted to az/el; flag = Impulse5 triangle.

```
in-train clicks: 135200 (files: 251006_001_0001.WAV, 251006_001_0002.WAV, 251006_001_0003.WAV)
flagged (closure>1 or infeasible): 60.8%

1. consecutive interval (ms), n=135197
     0.40-0.50        8093    6.0%
     0.50-0.75       20944   15.5%
     0.75-1.00       15223   11.3%
     1.00-1.50       21763   16.1%
     1.50-2.00       19950   14.8%
     2.00-2.50       20108   14.9%
     2.50-3.00        3586    2.7%
     3.00-5.00        7039    5.2%
     5.00-7.50        3659    2.7%
     7.50-10.00       2764    2.0%
    10.00-15.00       3617    2.7%
    15.00-20.00       2499    1.8%
    20.00-30.00       2465    1.8%
    30.00-50.00       1887    1.4%
    50.00-75.00        562    0.4%
    75.00-100.00       193    0.1%
   100.00-200.00       195    0.1%
   200.00-500.00       209    0.2%
   500.00-9999.00      441    0.3%

2. all-pair lag histogram, counts per 0.1 ms (row = 1 ms)
   0.0- 1.0 ms:      0      0      0      0   8093   9362   8265   6716   6363   6410
   1.0- 2.0 ms:   6969   6705   6556   6960   6990   7251   5949   6705  10515  15810
   2.0- 3.0 ms:  16323  11034   7480   6654   6693   7422   6565   6284   5467   5409
   3.0- 4.0 ms:   6271   6123   6251   5456   5301   6283   7557   8352   8906  10737
   4.0- 5.0 ms:   9906   9219   8437   7596   6789   6670   6025   5763   5523   5349
   5.0- 6.0 ms:   5288   5531   5749   5285   5901   6900   6878   7082   7546   8700
   6.0- 7.0 ms:   8168   7788   7578   7398   6793   6629   5941   5949   5630   5384
   7.0- 8.0 ms:   5691   5704   5741   5993   6195   6464   6322   6736   7205   7635
   8.0- 9.0 ms:   7289   7019   6939   6882   6905   6851   6106   5869   5728   5794
   9.0-10.0 ms:   6066   6097   6192   6147   6143   6283   6221   6544   7006   7053
  10.0-11.0 ms:   6869   6704   6661   6429   6847   7132   6117   6221   6054   5997
  11.0-12.0 ms:   6837   6326   6124   6501   6280   6612   6133   6597   7042   6774
  12.0-13.0 ms:   6637   6591   6708   6678   6558   7076   6563   6276   6690   6463
  13.0-14.0 ms:   6805   6355   6328   6315   6508   6738   6425   6946   6799   6591
  14.0-15.0 ms:   6637   6356   6618   6602   6736   6765   6747   6704   6453   6536

3. consecutive pairs by lag: amplitude ratio later/earlier, direction change (SRP, GCC; both unflagged SRP)
   0.40-0.75    n= 29037  weaker 54.3% ratio med 0.96 p10 0.43 | SRP med  46.4 <2deg  0.8% | GCC med  46.1 | unflagged n= 3760 med  54.5 <2deg  1.8%
   0.75-1.50    n= 36986  weaker 55.4% ratio med 0.93 p10 0.22 | SRP med  43.5 <2deg  1.5% | GCC med  42.8 | unflagged n= 5126 med  42.6 <2deg  4.6%
   1.50-1.80    n=  8911  weaker 59.1% ratio med 0.93 p10 0.21 | SRP med  42.8 <2deg  6.8% | GCC med  42.3 | unflagged n= 1654 med  50.2 <2deg 18.7%
   1.80-2.20    n= 24734  weaker 51.2% ratio med 1.00 p10 0.74 | SRP med  10.8 <2deg 35.4% | GCC med  13.0 | unflagged n= 7878 med   0.9 <2deg 73.6%
   2.20-3.50    n= 12239  weaker 50.6% ratio med 1.00 p10 0.48 | SRP med  28.4 <2deg 22.0% | GCC med  27.6 | unflagged n= 3261 med   1.3 <2deg 57.7%
   3.50-4.50    n=  3815  weaker 44.7% ratio med 1.03 p10 0.63 | SRP med  35.2 <2deg 13.7% | GCC med  34.6 | unflagged n=  737 med   6.3 <2deg 41.0%
   4.50-10.00   n=  7407  weaker 39.7% ratio med 1.05 p10 0.68 | SRP med  36.1 <2deg 14.5% | GCC med  37.5 | unflagged n= 1624 med   5.3 <2deg 45.4%
  10.00-30.00   n=  8581  weaker 33.4% ratio med 1.22 p10 0.59 | SRP med  42.0 <2deg 10.9% | GCC med  42.7 | unflagged n= 2123 med  33.1 <2deg 32.9%
  30.00-9999.00 n=  3487  weaker 36.5% ratio med 1.19 p10 0.58 | SRP med  44.1 <2deg  3.1% | GCC med  45.4 | unflagged n=  633 med  47.3 <2deg 11.8%
  random pairs SRP: median 45.5 deg, <5deg 1.3%

4. follower (0.4-1.6 ms after previous) vs predecessor, both unflagged
  n=9549 |daz|<3: 30.1% (random 7.4%); same az & |del|>5: 24.6%
  del p10/p25/p50/p75/p90: -62.9/-34.6/4.1/41.1/60.8

5. populations and flag rates
  isolated (gap>=1.6ms)  n= 65685 (48.6%) flagged 52.0% median ppeak 839
  follower               n= 69515 (51.4%) flagged 69.1% median ppeak 862
  amplitude Q1 weakest   flagged 74.5%
  amplitude Q2           flagged 69.6%
  amplitude Q3           flagged 61.9%
  amplitude Q4 strongest flagged 37.1%
  clean (isolated + unflagged): 31525 (23.3%)
  consecutive clean clicks, direction change by lag:
       1.6-2.2      n=  7674 median   0.9 deg <2deg  71.4%
       2.2-5.0      n=  6501 median   1.6 deg <2deg  54.6%
       5.0-10.0     n=  4557 median   2.7 deg <2deg  45.6%
      10.0-30.0     n=  8105 median   2.3 deg <2deg  48.6%
      30.0-100.0    n=  3541 median  11.1 deg <2deg  34.2%
     100.0-9999.0   n=  1144 median  32.8 deg <2deg   7.0%
```
