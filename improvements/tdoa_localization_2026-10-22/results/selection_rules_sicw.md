# Selection rules on real data (sicw)

```
feature table: click_features_sicw_New_device.csv (135200 clicks); consecutive isolated clicks (no detection 0.4-1.6 ms before either)
rule                                              ≠ rank 1 | 1.6-2.2 ms       | 2.2-10 ms        | 10-30 ms         |     in
                                                           | median / <2°     | median / <2°     | median / <2°     |  frame
rank 1 (SRP maximum)                                  0.0% |  11.2° /  33.8% |  23.6° /  24.6% |  31.3° /  20.2% |  61.3%
earliest arrival                                     51.0% |  21.0° /  25.6% |  27.1° /  20.0% |  31.6° /  17.2% |  62.9%
earliest, arriving <= 20 smp before the detection     39.2% |  18.4° /  28.2% |  26.9° /  21.5% |  33.7° /  17.4% |  59.5%
earliest, arriving <= 40 smp before the detection     44.9% |  19.5° /  27.2% |  26.7° /  20.6% |  32.9° /  17.2% |  60.6%
earliest, rscore >= 0.85 x best                      34.4% |  15.1° /  30.4% |  24.3° /  22.8% |  30.2° /  18.7% |  62.8%
earliest, rscore >= 0.95 x best                      19.9% |  12.2° /  32.6% |  23.0° /  24.0% |  30.3° /  19.5% |  62.5%
largest energy removed (efrac)                       33.0% |  15.1° /  31.2% |  26.7° /  22.6% |  34.2° /  18.3% |  57.5%
track lag rule                                       30.4% |  11.5° /  30.0% |  19.2° /  23.2% |  25.3° /  20.2% |  67.7%
pairs: 23003 (1.6-2.2 ms), 16082 (2.2-10 ms), 4667 (10-30 ms)
```
