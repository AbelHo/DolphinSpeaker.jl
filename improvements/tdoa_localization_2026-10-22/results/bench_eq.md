# Round-5 benchmark: hydrophone equalization (simulation)

```
15 scenarios x 4 s, 2 dolphins, calibrated realism (A0 6000, scatter 0.7, reverb 0.1); SIC -110:50, longest track + smart inherit, separation n = 2 (tuned)

version             | >=2 cand closure<=1     ncc | acc3° all acc3° dir.  med err |  sep F1 | phase err   |delay|
none                |    71.3%     44.2%    0.40 |    54.2%      58.1%     1.1° |   0.806 |      NaN      NaN
aligned             |    38.5%     83.6%    0.61 |     3.0%       4.1%    33.8° |   0.961 |     1.02      2.5
blind               |    49.3%     79.4%    0.65 |     0.4%       0.7%    30.0° |   0.916 |     0.91      NaN
blind + true delays |    30.8%     84.4%    0.66 |    34.0%      39.8%     6.4° |   0.943 |     0.91      NaN
oracle              |    23.0%     92.4%    0.81 |    91.0%      98.0%     0.3° |   0.971 |     0.00      0.0

aligned             - none: accuracy -60.3 points (median; better in 1 of 15), separation F1 +0.023 (median); share of the oracle's gain recovered -366% (scenarios where the oracle gains >= 5 points)
blind               - none: accuracy -60.8 points (median; better in 0 of 15), separation F1 +0.022 (median); share of the oracle's gain recovered -370% (scenarios where the oracle gains >= 5 points)
blind + true delays - none: accuracy -6.1 points (median; better in 6 of 15), separation F1 +0.023 (median); share of the oracle's gain recovered -287% (scenarios where the oracle gains >= 5 points)
```
