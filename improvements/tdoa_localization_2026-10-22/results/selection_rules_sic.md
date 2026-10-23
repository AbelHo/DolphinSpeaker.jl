# Selection rules on real data (sic)

```
feature table: click_features_sic_New_device.csv (135200 clicks); consecutive isolated clicks (no detection 0.4-1.6 ms before either)
rule                                              ≠ rank 1 | 1.6-2.2 ms       | 2.2-10 ms        | 10-30 ms         |     in
                                                           | median / <2°     | median / <2°     | median / <2°     |  frame
rank 1 (SRP maximum)                                  0.0% |  10.4° /  35.6% |  24.3° /  25.1% |  33.6° /  19.5% |  57.4%
earliest arrival                                     49.8% |  20.7° /  26.7% |  27.3° /  20.2% |  34.8° /  16.2% |  59.4%
earliest, arriving <= 20 smp before the detection     40.6% |  19.0° /  28.7% |  28.0° /  21.4% |  36.2° /  16.8% |  56.2%
earliest, arriving <= 40 smp before the detection     47.8% |  20.5° /  27.1% |  27.9° /  20.4% |  35.6° /  16.3% |  58.0%
earliest, rscore >= 0.85 x best                      34.0% |  14.7° /  31.3% |  24.8° /  22.8% |  33.5° /  18.1% |  59.1%
earliest, rscore >= 0.95 x best                      20.3% |  11.9° /  33.7% |  24.1° /  24.1% |  32.9° /  19.2% |  58.8%
largest energy removed (efrac)                       32.6% |  15.2° /  32.0% |  27.3° /  23.2% |  35.5° /  17.5% |  53.9%
track lag rule                                       30.5% |  11.6° /  31.1% |  19.6° /  23.2% |  28.6° /  19.2% |  63.9%
pairs: 23003 (1.6-2.2 ms), 16082 (2.2-10 ms), 4667 (10-30 ms)
```
