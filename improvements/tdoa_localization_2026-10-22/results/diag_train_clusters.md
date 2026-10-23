# Direction clusters per click train

Source: `/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device`. Clean = gap to previous detection >= 1.6 ms and Impulse5 not flagged; directions from Impulse6 (SRP).

```
clean clicks: 31525 of 135200

2D histogram of clean SRP directions (% of clean, 10 deg bins; rows el, cols az)
   el\az  -90  -80  -70  -60  -50  -40  -30  -20  -10    0   10   20   30   40   50   60   70   80
      80    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .
      70    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .
      60    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .
      50    .    .    .    .    .    .  0.4  0.5    .  0.4    .  0.9    .    .    .    .    .    .
      40    .    .    .    .    .    .    .    .  0.4  1.1  0.5    .    .    .    .    .    .    .
      30    .    .    .    .    .    .    .    .  0.3  1.0  0.4  1.5    .    .  0.3    .    .    .
      20    .    .    .    .    .    .  0.4    .  0.3  0.3  0.6  1.9  1.0  0.5    .  0.3    .    .
      10    .    .    .    .    .    .    .  0.4  1.2  1.4  1.7  0.9  0.7    .  0.5  0.6  2.3    .
       0    .    .    .    .    .    .    .  0.8  1.3  3.0  2.0  0.9  0.5  0.5  0.7  0.7    .    .
     -10    .    .    .    .    .  0.6  0.5  0.7  3.9  5.1  1.5  1.1  0.8  0.5  0.4  0.5    .    .
     -20    .    .    .    .    .    .  0.5  1.9  5.2  5.5  1.7  1.4  0.4  0.4  0.4    .    .    .
     -30    .    .    .    .    .    .    .  1.4  2.3  2.5  3.0  1.4    .  0.6  1.2  1.0  0.3    .
     -40    .    .    .    .    .    .    .    .  0.4  0.9  0.5    .  0.3    .    .    .    .    .
     -50    .    .    .    .    .    .    .    .  0.6    .  0.3  0.3  0.4  0.4  0.3    .    .    .
     -60    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .
     -70    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .
     -80    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .
     -90    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .    .

trains analysed (>= 20 clean clicks): 215
  0 clusters >= 10%: 0 trains
  1 clusters >= 10%: 56 trains
  2 clusters >= 10%: 80 trains
  3 clusters >= 10%: 63 trains
  4 clusters >= 10%: 14 trains
  5 clusters >= 10%: 1 trains
  6 clusters >= 10%: 1 trains
  share of largest cluster: median 0.50; second: median 0.17

second cluster (n trains=159): daz p25/p50/p75 -10.8/0.4/13.2, del -18.5/2.3/34.3, |daz|<5: 31%
lag from the latest cluster-1 detection to each clean cluster-2 click (random times in the same train):
    0.0-1.6   ms:   0.0%  (random  11.0%)
    1.6-2.5   ms:   8.0%  (random   3.8%)
    2.5-5.0   ms:   6.2%  (random   5.6%)
    5.0-10.0  ms:   7.3%  (random   8.0%)
   10.0-30.0  ms:  14.1%  (random  13.2%)
   30.0-Inf   ms:  53.9%  (random  50.3%)
ICI (ms) cluster1 p10/p25/p50/p75/p90: 1.98/2.10/5.88/18.90/43.91
ICI (ms) cluster2 p10/p25/p50/p75/p90: 2.06/4.05/13.06/26.90/63.76
```
