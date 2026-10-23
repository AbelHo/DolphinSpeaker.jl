# Dolphin separation on New_device (tuned parameters)

```
New_device, 135200 detections; directions: longest-track rule on SIC -110:50 candidates; separate_sources parameters: tuned (sigma_meas_deg = 8.0, gate_deg = 25.0, force_gate_deg = 40.0, t_drop = 1.0, follow_gate_deg = 12.0, merge_deg = 8.0)

=== n = 2: unassigned 11.6%; roles: leader 45.2%, own 5.9%, repeat 6.2%, echo 31.1%, unassigned 11.6%
  label  detections   az med   el med    ICI p25    ICI med    ICI p75
  1          67287     3.1°    -1.8°      2.0 ms      3.3 ms      8.3 ms
  2          52220     2.6°    -0.5°      2.1 ms      4.0 ms     11.0 ms
  consecutive genuine clicks 2-30 ms apart: within a dolphin median change 5.5° (<2°: 33.7%), ICI change 20% | unseparated 11.0° (<2°: 27.7%), ICI change 15%
  inside the frame: 72.7% (own directions 68.0%)

=== n estimated: unassigned 4.5%; roles: leader 47.6%, own 15.8%, repeat 7.3%, echo 24.8%, unassigned 4.5%
  confirmed sources per file: 0001 370, 0002 582, 0003 182
  simultaneously active sources per 1 s with clicks: 0: 3%, 1: 22%, 2: 23%, 3: 18%, 4: 13%
  label  detections   az med   el med    ICI p25    ICI med    ICI p75
  3038        2620     0.0°     4.2°      1.4 ms      3.4 ms      7.7 ms
  491         2521     1.8°     4.6°      1.3 ms      2.0 ms      3.7 ms
  1080        2235     3.2°     0.5°      1.6 ms      2.7 ms      5.9 ms
  2389        2187     8.2°   -11.7°      2.0 ms      2.1 ms      3.3 ms
  726         2157     3.4°    -8.6°      1.5 ms      2.6 ms      5.2 ms
  576         2113     5.1°     3.5°      1.8 ms      2.2 ms      3.7 ms
  consecutive genuine clicks 2-30 ms apart: within a dolphin median change 3.8° (<2°: 36.8%), ICI change 20% | unseparated 14.7° (<2°: 24.6%), ICI change 11%
  inside the frame: 72.7% (own directions 68.0%)
```
