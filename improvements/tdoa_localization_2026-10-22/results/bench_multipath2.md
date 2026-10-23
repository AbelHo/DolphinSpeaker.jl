# Round-2 multipath benchmark: candidate generators

```
trials 3000 (a third without echo), templates 300 real clean clicks, real noise, SNR 19 dB, window -50:50 (sic_wide: -110:50)
recall = direct path within 2 deg among the accepted candidates (echo trials); echo found = an accepted extra candidate within 3 deg of an echo;
false = accepted extra candidates within 3 deg of no true arrival (per click); rank1 / earliest = selection by rank 1 / by earliest beam arrival (echo trials)

method  rel>=  supp>= |  recall   rank1 earliest | echo found  false/click false/click
                      |  (echo)                  |               (no echo)     (echo)
plain   0.0    -      |    67.9    55.7     48.7 |       18.6         2.00       1.67
plain   0.3    -      |    67.9    55.7     48.7 |       18.6         2.00       1.67
plain   0.5    -      |    67.9    55.7     49.3 |       18.5         1.60       1.60
plain   0.7    -      |    67.1    55.7     58.2 |       15.6         0.09       0.78
plain   0.9    -      |    60.5    55.7     58.4 |        7.0         0.00       0.19
plain   2.0    -      |    55.7    55.7     55.7 |        0.0         0.00       0.00
plain   0.0    0.3    |    67.3    55.7     53.4 |       17.3         0.60       1.13
plain   0.3    0.3    |    67.3    55.7     53.4 |       17.3         0.60       1.13
plain   0.5    0.3    |    67.3    55.7     53.4 |       17.3         0.60       1.13
plain   0.7    0.3    |    66.8    55.7     58.5 |       15.4         0.09       0.73
plain   0.9    0.3    |    60.4    55.7     58.3 |        7.0         0.00       0.19
plain   2.0    0.3    |    55.7    55.7     55.7 |        0.0         0.00       0.00
plain   0.0    0.5    |    65.7    55.7     59.0 |       14.3         0.08       0.62
plain   0.3    0.5    |    65.7    55.7     59.0 |       14.3         0.08       0.62
plain   0.5    0.5    |    65.7    55.7     59.0 |       14.3         0.08       0.62
plain   0.7    0.5    |    65.4    55.7     59.2 |       13.9         0.06       0.57
plain   0.9    0.5    |    60.3    55.7     58.2 |        7.0         0.00       0.18
plain   2.0    0.5    |    55.7    55.7     55.7 |        0.0         0.00       0.00
clean   0.0    -      |    64.2    55.7     41.4 |       13.7         2.00       1.76
clean   0.3    -      |    64.2    55.7     45.5 |       13.6         1.13       1.54
clean   0.5    -      |    62.7    55.7     56.3 |       11.0         0.11       0.79
clean   0.7    -      |    59.9    55.7     58.6 |        5.8         0.01       0.30
clean   0.9    -      |    57.2    55.7     57.2 |        1.9         0.00       0.06
clean   2.0    -      |    55.7    55.7     55.7 |        0.0         0.00       0.00
sic     0.0    -      |    75.0    55.7     58.5 |       39.3         2.00       1.37
sic     0.3    -      |    75.0    55.7     59.1 |       39.3         1.89       1.35
sic     0.5    -      |    74.7    55.7     66.9 |       38.4         0.32       0.87
sic     0.7    -      |    73.9    55.7     71.6 |       32.4         0.02       0.40
sic     0.9    -      |    70.7    55.7     70.3 |       23.7         0.00       0.15
sic     2.0    -      |    55.8    55.7     55.8 |        0.1         0.00       0.00
sic_wide 0.0    -      |    86.5    58.6     66.3 |       38.6         2.00       1.29
sic_wide 0.3    -      |    86.5    58.6     66.9 |       38.6         1.79       1.26
sic_wide 0.5    -      |    86.4    58.6     78.0 |       37.8         0.24       0.79
sic_wide 0.7    -      |    85.7    58.6     83.4 |       33.1         0.01       0.34
sic_wide 0.9    -      |    81.0    58.6     80.6 |       24.7         0.00       0.12
sic_wide 2.0    -      |    58.9    58.6     58.9 |        0.2         0.00       0.00
```
