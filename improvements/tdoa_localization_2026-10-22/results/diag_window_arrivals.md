# Arrivals inside the TDOA window

Source: `/media/spin/anas2/data/calf/upload/results/tmp_20261007/Two-way_2males_06_10_2025_multitdoa/New_device`, audio `/media/spin/anas2/data/calf/upload/Two-way/Two-way_2males_06_10_2025/New_device/251006_001_0001.WAV`. Snippets are high-passed like the pipeline (Butterworth, `impulsive_band_pass`). Ref channel = 3rd hydrophone column.

```
A. sampled clicks of 251006_001_0001.WAV: 3764 (every 15th)
envelope peaks (>= 30% of channel max) inside -50:50, ref channel:
  0:   1.9%
  1:  13.1%
  2:  18.8%
  3:  19.6%
  4+:  46.5%
>= 2 peaks on some channel: 92.9%
preceded by an earlier peak (8-150 samples): 83.6%; lag p25/p50/p75 78/119/132 samples; height med 0.48
  1 peak on every channel      n=  234 flagged  20.9%
  >= 2 peaks on some channel   n= 3497 flagged  61.0%
  not preceded                 n=  619 flagged  22.1%
  preceded                     n= 3145 flagged  65.7%
  isolated                     n= 1787 >=2 peaks  94.0% preceded  76.0% flagged  47.7%
  follower                     n= 1977 >=2 peaks  92.0% preceded  90.3% flagged  68.3%

B. 3000 consecutive clicks (3 blocks); earliest peak > 20 samples before detection: 97.1%
  earliest-peak vs detection direction: median 35.9 deg; closure<=1: detection 32%, earliest 22%
  consecutive isolated clicks, detection window -50:50      n=1724 median  14.9 deg <2deg  27.4%
  consecutive isolated clicks, earliest-peak window -15:40  n=1724 median  37.5 deg <2deg   5.1%
```
