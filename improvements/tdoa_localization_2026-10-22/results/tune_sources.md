# Tuning separate_sources on the realistic simulation

```
grid of 288 parameter sets, tuned on 27 two-dolphin scenarios, checked on 14 three-dolphin scenarios (not used for tuning)
best for n given: (sigma_meas_deg = 8.0, gate_deg = 25.0, force_gate_deg = 40.0, t_drop = 1.0, follow_gate_deg = 12.0, merge_deg = 4.0)
best for n estimated: (sigma_meas_deg = 8.0, gate_deg = 25.0, force_gate_deg = 40.0, t_drop = 1.0, follow_gate_deg = 12.0, merge_deg = 8.0)

                                       F1 precision  recall  acc3° unassigned%  n found
2 dolphins (tuning set)
  defaults, n given                 0.777     0.966   0.691   44.4        16.3      1.7
  defaults, n estimated             0.787     0.968   0.701   44.5         4.9      5.2
  tuned, n given                    0.839     0.949   0.777   44.3         6.3      1.8
  tuned (for estimated), n estimated  0.829     0.956   0.759   44.2         2.1      4.0
  k-means baseline (k = n)          0.851
3 dolphins (held out)
  defaults, n given                 0.769     0.928   0.697   54.7        17.6      2.8
  defaults, n estimated             0.797     0.930   0.731   54.3         3.3      7.8
  tuned, n given                    0.828     0.913   0.778   54.5         1.9      2.8
  tuned (for estimated), n estimated  0.832     0.919   0.778   54.5         1.3      4.9
  k-means baseline (k = n)          0.801

top 10 parameter sets (n given, 2 dolphins):
10×13 DataFrame
 Row │ sigma_meas  gate     force_gate  t_drop   follow_gate  merge    f1_given  f1_est    recall_given  precision_given  acc_given  unassigned_given  n_est
     │ Float64     Float64  Float64     Float64  Float64      Float64  Float64   Float64   Float64       Float64          Float64    Float64           Float64
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │        8.0     25.0        40.0      1.0         12.0      4.0  0.838764  0.827148      0.777431         0.949248    44.2973           6.29212  4.07407
   2 │        8.0     25.0        40.0      1.0         12.0      8.0  0.838764  0.828509      0.777431         0.949248    44.2973           6.29212  4.03704
   3 │        8.0     25.0        40.0      3.0         12.0      4.0  0.837509  0.821979      0.775899         0.948272    44.282            5.91815  3.48148
   4 │        8.0     25.0        40.0      3.0         12.0      8.0  0.837509  0.823565      0.775899         0.948272    44.282            5.91815  3.37037
   5 │        8.0     25.0        40.0      1.0          5.0      4.0  0.835714  0.826451      0.77408          0.948845    44.5377           6.64742  4.07407
   6 │        8.0     25.0        40.0      1.0          5.0      8.0  0.835714  0.827812      0.77408          0.948845    44.5377           6.64742  4.03704
   7 │       12.0     25.0        40.0      1.0         12.0      4.0  0.835225  0.827306      0.774117         0.946505    44.2717           6.50167  4.03704
   8 │       12.0     25.0        40.0      1.0         12.0      8.0  0.835225  0.827306      0.774117         0.946505    44.2717           6.50167  4.03704
   9 │        8.0     25.0        40.0      3.0          5.0      4.0  0.834722  0.820462      0.77301          0.947772    44.5469           6.27834  3.51852
  10 │        8.0     25.0        40.0      3.0          5.0      8.0  0.834722  0.82405       0.77301          0.947772    44.5469           6.27834  3.33333
```
