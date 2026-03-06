# Kassab et al. 1993 — Real Data Extraction

Source: "Morphometry of pig coronary arterial trees" (Am J Physiol Heart 265:H350-H365)
PDF: `/Users/daleblack/Documents/vessel tree growing papers/kassab-et-al-1993-morphometry-of-pig-coronary-arterial-trees.pdf`

## Key Concepts

### Segments vs Elements
- A **segment** is the vessel between two consecutive bifurcation nodes
- An **element** is a series of segments of the SAME Strahler order connected end-to-end
- Element diameter = average of constituent segment diameters
- Element length = sum of constituent segment lengths
- The S/E ratio (Table 5) gives the average number of segments per element
- **The connectivity matrix operates on ELEMENTS, not segments**

### Diameter-Defined Strahler Ordering
Kassab does NOT use pure topological Strahler ordering. The method (detailed in Jiang 1994) is:
1. Start with conventional Strahler ordering
2. Compute mean D_n and SD_n per order
3. Define non-overlapping diameter ranges via:
   - Lower bound: D'_1(n) = [(D_{n-1} + SD_{n-1}) + (D_n - SD_n)] / 2
   - Upper bound: D'_2(n) = [(D_n + SD_n) + (D_{n+1} - SD_{n+1})] / 2
4. Reassign orders based on diameter ranges
5. Iterate until convergence (<1% of vessels change order)

See `jiang_1994_methodology.md` for the full algorithm.

---

## Table 1: RCA Segments and Elements

### RCA Segments (Diameter and Length in um)

| Order | N_seg    | D_seg(um) | SD_D   | L_seg(um) | SD_L   |
|-------|----------|-----------|--------|-----------|--------|
| 1     | 854,880  | 10.8      | 2.85   | 96        | 78     |
| 2     | 163,233  | 21.5      | 3.48   | 92        | 74     |
| 3     | 43,858   | 39.1      | 5.64   | 100       | 83     |
| 4     | 12,285   | 69.2      | 10.5   | 123       | 111    |
| 5     | 5,449    | 124       | 20.2   | 176       | 166    |
| 6     | 3,279    | 228       | 40.4   | 311       | 299    |
| 7     | 1,563    | 408       | 68.8   | 503       | 440    |
| 8     | 550      | 703       | 118    | 810       | 720    |
| 9     | 135      | 1211      | 210    | 1309      | 1260   |
| 10    | 36       | 2167      | 290    | 1997      | 1750   |
| 11    | 6        | 3480      | 440    | 2610      | 1870   |

### RCA Elements (Diameter and Length in um)

| Order | N_elem   | D_elem(um) | SD_D   | L_elem(um) | SD_L   |
|-------|----------|------------|--------|------------|--------|
| 1     | 393,294  | 11.1       | 2.88   | 231        | 211    |
| 2     | 80,968   | 21.7       | 3.46   | 196        | 184    |
| 3     | 20,898   | 39.3       | 5.75   | 220        | 218    |
| 4     | 5,753    | 69.7       | 10.3   | 277        | 289    |
| 5     | 2,469    | 125        | 19.9   | 411        | 408    |
| 6     | 1,275    | 231        | 40.0   | 872        | 1040   |
| 7     | 518      | 413        | 72.8   | 1629       | 1950   |
| 8     | 175      | 714        | 118    | 2770       | 3560   |
| 9     | 46       | 1210       | 233    | 5620       | 7660   |
| 10    | 13       | 2233       | 287    | 8540       | 8020   |
| 11    | 3        | 3637       | 407    | 10070      | 10100  |

---

## Table 2: LAD Segments and Elements

### LAD Segments

| Order | N_seg    | D_seg(um) | SD_D   | L_seg(um) | SD_L   |
|-------|----------|-----------|--------|-----------|--------|
| 1     | 823,990  | 10.6      | 2.68   | 97        | 79     |
| 2     | 155,063  | 21.5      | 3.66   | 93        | 76     |
| 3     | 42,710   | 39.3      | 5.93   | 104       | 88     |
| 4     | 10,820   | 71.2      | 11.1   | 129       | 121    |
| 5     | 5,108    | 129       | 21.3   | 194       | 182    |
| 6     | 3,302    | 228       | 36.6   | 310       | 305    |
| 7     | 1,432    | 395       | 59.5   | 471       | 434    |
| 8     | 431      | 676       | 107    | 801       | 682    |
| 9     | 104      | 1200      | 207    | 1460      | 1480   |
| 10    | 18       | 2040      | 354    | 2820      | 3430   |
| 11    | 3        | 3307      |        | 3320      |        |

### LAD Elements

| Order | N_elem   | D_elem(um) | SD_D   | L_elem(um) | SD_L   |
|-------|----------|------------|--------|------------|--------|
| 1     | 368,554  | 10.8       | 2.72   | 230        | 210    |
| 2     | 77,267   | 21.6       | 3.64   | 198        | 189    |
| 3     | 20,370   | 39.7       | 5.98   | 230        | 235    |
| 4     | 5,091    | 71.8       | 11.2   | 291        | 313    |
| 5     | 2,255    | 130        | 21.0   | 463        | 462    |
| 6     | 1,198    | 233        | 38.3   | 906        | 1130   |
| 7     | 470      | 401        | 62.9   | 1530       | 1890   |
| 8     | 142      | 688        | 107    | 2700       | 3200   |
| 9     | 35       | 1210       | 233    | 6330       | 8640   |
| 10    | 9        | 2207       | 298    | 7580       | 7350   |
| 11    | 2        | 3307       |        | 10000      |        |

---

## Table 3: LCX Segments and Elements

### LCX Segments

| Order | N_seg    | D_seg(um) | SD_D   | L_seg(um) | SD_L   |
|-------|----------|-----------|--------|-----------|--------|
| 1     | 518,093  | 10.3      | 2.51   | 93        | 76     |
| 2     | 107,853  | 20.1      | 3.52   | 89        | 73     |
| 3     | 27,568   | 37.1      | 5.63   | 94        | 80     |
| 4     | 7,368    | 66.3      | 9.97   | 113       | 102    |
| 5     | 4,043    | 120       | 19.0   | 178       | 167    |
| 6     | 2,403    | 215       | 34.9   | 290       | 272    |
| 7     | 953      | 371       | 59.6   | 433       | 380    |
| 8     | 310      | 636       | 103    | 747       | 710    |
| 9     | 54       | 1113      | 184    | 1130      | 960    |
| 10    | 9        | 1990      | 300    | 2020      | 1750   |

### LCX Elements

| Order | N_elem   | D_elem(um) | SD_D   | L_elem(um) | SD_L   |
|-------|----------|------------|--------|------------|--------|
| 1     | 244,025  | 10.5       | 2.55   | 206        | 186    |
| 2     | 54,463   | 20.3       | 3.48   | 185        | 177    |
| 3     | 13,368   | 37.3       | 5.68   | 202        | 207    |
| 4     | 3,420    | 66.7       | 10.2   | 253        | 275    |
| 5     | 1,688    | 121        | 19.1   | 459        | 453    |
| 6     | 888      | 218        | 35.6   | 820        | 940    |
| 7     | 321      | 379        | 61.1   | 1380       | 1500   |
| 8     | 101      | 654        | 104    | 2530       | 3060   |
| 9     | 21       | 1137       | 218    | 4090       | 4380   |
| 10    | 5        | 2020       | 300    | 5670       | 5070   |

---

## Table 5: S/E Ratios (Segments per Element)

| Order | RCA   | LAD   | LCX   |
|-------|-------|-------|-------|
| 1     | 2.17  | 2.24  | 2.12  |
| 2     | 2.02  | 2.01  | 1.98  |
| 3     | 2.10  | 2.10  | 2.06  |
| 4     | 2.14  | 2.13  | 2.15  |
| 5     | 2.21  | 2.26  | 2.40  |
| 6     | 2.57  | 2.76  | 2.71  |
| 7     | 3.02  | 3.05  | 2.97  |
| 8     | 3.14  | 3.04  | 3.07  |
| 9     | 2.93  | 2.97  | 2.57  |
| 10    | 2.77  | 2.00  | 1.80  |
| 11    | 2.00  | 1.50  | —     |

---

## Table 6: Connectivity Matrix — RCA (REAL DATA)

CM[daughter_order, parent_order] = mean number of elements of daughter_order
that spring directly from elements of parent_order, divided by total elements of parent_order.

Format: value +/- SE. Unlisted entries are 0.

```
Parent order →  1        2        3        4        5        6        7        8        9        10       11
Daughter ↓
0:           2.75±.082 0.674±.067 0.151±.045 0.040±.028
1:           0.131±.027 2.13±.070 0.802±.101 0.300±.104 0.008±.008 0.004±.003
2:                     0.080±.025 2.15±.085 0.700±.109 0.159±.037 0.020±.006 0.037±.013
3:                                0.070±.028 2.12±.116 0.688±.067 0.314±.027 0.244±.033 0.143±.033 0.059±.041
4:                                           0.160±.052 2.13±.090 0.444±.031 0.640±.057 0.468±.074 0.324±.145
5:                                                      0.344±.035 2.43±.046 1.43±.074  1.51±.157  0.853±.189 0.727±.304
6:                                                                 0.167±.018 2.14±.060 1.37±.106  1.68±.334  1.80±.611  1
7:                                                                            0.130±.019 2.32±.092 1.23±.246  2.80±.854  11
8:                                                                                       0.099±.031 2.62±.227 2.00±.907  8
9:                                                                                                  0.059±.041 2.40±.562 5
10:                                                                                                            0.400±.221 6
```

### RCA CM as code-ready matrix (12x12, 0-indexed orders → 1-indexed Julia):
```julia
CM = zeros(Float64, 12, 12)
# Column = parent order+1, Row = daughter order+1
# Parent order 1 (col 2)
CM[1,2]=2.75; CM[2,2]=0.131
# Parent order 2 (col 3)
CM[1,3]=0.674; CM[2,3]=2.13; CM[3,3]=0.080
# Parent order 3 (col 4)
CM[1,4]=0.151; CM[2,4]=0.802; CM[3,4]=2.15; CM[4,4]=0.070
# Parent order 4 (col 5)
CM[1,5]=0.040; CM[2,5]=0.300; CM[3,5]=0.700; CM[4,5]=2.12; CM[5,5]=0.160
# Parent order 5 (col 6)
CM[2,6]=0.008; CM[3,6]=0.159; CM[4,6]=0.688; CM[5,6]=2.13; CM[6,6]=0.344
# Parent order 6 (col 7)
CM[2,7]=0.004; CM[3,7]=0.020; CM[4,7]=0.314; CM[5,7]=0.444; CM[6,7]=2.43; CM[7,7]=0.167
# Parent order 7 (col 8)
CM[3,8]=0.037; CM[4,8]=0.244; CM[5,8]=0.640; CM[6,8]=1.43; CM[7,8]=2.14; CM[8,8]=0.130
# Parent order 8 (col 9)
CM[4,9]=0.143; CM[5,9]=0.468; CM[6,9]=1.51; CM[7,9]=1.37; CM[8,9]=2.32; CM[9,9]=0.099
# Parent order 9 (col 10)
CM[4,10]=0.059; CM[5,10]=0.324; CM[6,10]=0.853; CM[7,10]=1.68; CM[8,10]=1.23; CM[9,10]=2.62; CM[10,10]=0.059
# NOTE: CM[10,10] should be read carefully — the paper shows 0.059±0.041 for row 9 col 10
# Actually row 9 (order 9→daughter) col 10 (parent order 9) = 0.059
# Row 10 (order 10→daughter) col 10 (parent order 10) = 0.400
# Parent order 10 (col 11)
CM[5,11]=0.727; CM[6,11]=1.80; CM[7,11]=2.80; CM[8,11]=2.00; CM[9,11]=2.40; CM[10,11]=0.400; CM[11,11]=0.400
# NOTE: Need to recheck col 11 entries — some may be from parent order 11
# Parent order 11 (col 12)
CM[6,12]=1; CM[7,12]=11; CM[8,12]=8; CM[9,12]=5; CM[10,12]=6; CM[11,12]=0
```

**IMPORTANT**: The order-11 column values (1, 11, 8, 5, 6) represent total counts from typically 1-3 trees, NOT ratios. They should be normalized by N_parent(11) to get per-parent ratios. For RCA, N_elem(11)=3, so divide by 3.

---

## Table 7: Connectivity Matrix — LAD (REAL DATA)

```
Parent order →  1        2        3        4        5        6        7        8        9        10       11
Daughter ↓
0:           3.18±.118 0.675±.080 0.148±.081
1:           0.144±.031 2.04±.070 0.630±.116 0.071±.071
2:                     0.094±.027 2.24±.102 1.50±.274  0.063±.025 0.094±.019 0.023±.010
3:                                0.074±.036 2.14±.204 0.381±.056 0.098±.018 0.120±.028 0.092±.029 0.030±.030
4:                                           0.143±.097 2.25±.097 0.425±.042 0.380±.047 0.428±.070 0.303±.102 0.167±.167
5:                                                      0.238±.035 2.50±.070 1.91±.111  1.58±.173  1.38±.281  0.667±.333
6:                                                                 0.155±.022 2.50±.080 1.58±.159  1.48±.323  1.17±.654
7:                                                                            0.116±.022 2.09±.106 1.30±.300  2.00±1.09  2
8:                                                                                       0.061±.024 2.50±.209 2.50±1.06  3
9:                                                                                                  0.121±.058 3.33±1.12 8
10:                                                                                                            0.100±.100 5
```

---

## Table 8: Connectivity Matrix — LCX (REAL DATA)

```
Parent order →  1        2        3        4        5        6        7        8        9        10
Daughter ↓
0:           3.18±.118 0.675±.080 0.148±.081
1:           0.144±.031 2.04±.070 0.630±.116 0.071±.071
2:                     0.094±.027 2.24±.102 1.50±.274  0.150±.058
3:                                0.074±.036 2.14±.204 0.150±.077 0.025±.018 0.011±.011
4:                                           0.143±.097 2.85±.148 0.385±.076 0.179±.042 0.109±.046
5:                                                      0.100±.043 2.51±.089 1.13±.115  0.956±.158 0.444±.242
6:                                                                 0.213±.042 2.33±.110 2.02±.252  1.78±.325  2
7:                                                                            0.168±.046 2.04±.132 2.00±.667  4
8:                                                                                       0.196±.067 4.00±.680 4
9:                                                                                                  0.111±.111 8
```

**NOTE**: LCX has only 10 orders (0-9 as rows, 1-10 as parent columns), not 12.

---

## Table 9: Total Number of Vessel Elements per Order (mean +/- SE across specimens)

### RCA
| Order | N_elements      |
|-------|-----------------|
| 1     | 393,294±158,657 |
| 2     | 80,968±34,753   |
| 3     | 20,898±8,935    |
| 4     | 5,753±2,446     |
| 5     | 2,469±1,059     |
| 6     | 1,275±553       |
| 7     | 518±203         |
| 8     | 175±52          |
| 9     | 46±16           |
| 10    | 13±5            |
| 11    | 3±1             |

### LAD
| Order | N_elements      |
|-------|-----------------|
| 1     | 368,554±221,134 |
| 2     | 77,267±44,768   |
| 3     | 20,370±11,838   |
| 4     | 5,091±2,832     |
| 5     | 2,255±1,393     |
| 6     | 1,198±780       |
| 7     | 470±278         |
| 8     | 142±66          |
| 9     | 35±13           |
| 10    | 9±3             |
| 11    | 2±1             |

### LCX
| Order | N_elements      |
|-------|-----------------|
| 1     | 244,025±121,403 |
| 2     | 54,463±25,937   |
| 3     | 13,368±5,937    |
| 4     | 3,420±1,377     |
| 5     | 1,688±736       |
| 6     | 888±360         |
| 7     | 321±145         |
| 8     | 101±30          |
| 9     | 21±7            |
| 10    | 5±2             |

---

## Table 4: Empirical Constants for Log-Linear Fits

Kassab fits: log10(Y) = a + b * order

### Diameter: log10(D_n) = a + b*n
| Artery | a (segments) | b (segments) | a (elements) | b (elements) |
|--------|-------------|-------------|-------------|-------------|
| RCA    | 0.790       | 0.243       | 0.802       | 0.243       |
| LAD    | 0.783       | 0.242       | 0.792       | 0.243       |
| LCX    | 0.773       | 0.241       | 0.782       | 0.242       |

### Length: log10(L_n) = a + b*n
| Artery | a (segments) | b (segments) | a (elements) | b (elements) |
|--------|-------------|-------------|-------------|-------------|
| RCA    | 1.786       | 0.131       | 2.131       | 0.162       |
| LAD    | 1.787       | 0.132       | 2.124       | 0.165       |
| LCX    | 1.775       | 0.128       | 2.090       | 0.155       |

### Total number: log10(N_n) = a + b*n
| Artery | a (segments) | b (segments) | a (elements) | b (elements) |
|--------|-------------|-------------|-------------|-------------|
| RCA    | 6.328       | -0.436      | 6.006       | -0.441      |
| LAD    | 6.322       | -0.439      | 5.970       | -0.445      |
| LCX    | 6.108       | -0.430      | 5.786       | -0.436      |

---

## Key Observations for Implementation

1. **Three separate CMs needed**: RCA, LAD, and LCX have meaningfully different connectivity matrices. A single "coronary" CM is wrong.

2. **Element-level CM**: The CM operates on ELEMENTS, not segments. Generated trees must group segments into elements before computing CM statistics.

3. **Segment diameters differ from element diameters**: Segment means are slightly lower because they include partial-element segments at boundaries.

4. **Real diameter SDs are MUCH smaller than currently used**: e.g., Order 1 RCA segment SD = 2.85um (not 4.0), Order 5 = 20.2um (not 35.0).

5. **Real lengths are SHORT**: Order 1 segment length = 96um (not 50um), but element length = 231um. The CCO+subdivision approach should match element lengths since each "segment" in our tree represents an element between bifurcations.

6. **Highest order varies by artery**: RCA/LAD have 11 orders, LCX has 10.

7. **Asymmetry**: Not explicitly in the 1993 paper tables. The Beta(2.5, 0.8) comes from later Kassab work. Need to verify.

8. **Order-0 (capillaries)**: Not in the tables above. Kassab excludes capillary network (not tree-like). Our order-0 should be precapillary arterioles.
