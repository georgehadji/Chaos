# Calibration tables: p3-abstention

record SHA-256: a8cfdbe12a56697fba17db12b26cafced06647fd05dfa0db79af9a1319b3680a

## Header

- preregistration header SHA-256 (recomputed under the rule of the SELF cell): f614a38449f611a42c56a0aa79bbc2d106dfaab679d1972a2b5349011df1a072
- p2 manifest SHA-256: cdef3cf423e16a0b54e81e96fd460f4cb32d7b788325ffbf4e0581f1db27acce
- P1 manifest SHA-256: ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d
- panel file SHA-256: 881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5
- p3 parameter table SHA-256: 69965ec24cee507d76924423745cf332ecbb6c003880e5fb4b1036a60fc6a5e1
- p2 parameter table SHA-256: f5278e0f5c499f06674cc23cbe93d7d0a23720f4e4f77ece6ab48c4659063099
- P1 parameter table SHA-256: a544130bc56dc66d8b1f7ec5d1e54197306f14dc19c7cd22fdbd94e4534b7282
- p2 preregistration SHA-256: a69bddb2a9045f8aac25ba7771fe7e5d90e70f9ff63203646e3684c8e19efff4
- P1 preregistration SHA-256: 62011ba5af514da4f1238aef9ae6cb35ca21093c4343b71dd696df360b9994bf
- module version: 0.1.2
- git commit: 02994688354738e11405605e8a0074481a7ccb80
- git working tree clean: true
- generator: xoshiro256++ seeded by SplitMix64 from the first 8 bytes, big-endian, of SHA-256 over chaoslab-analysis-1|<test_seed>|p3|<purpose>|<unit>
- quantile function: AS241 (PPND16)

## Table T1

### Relabelling C0 and the drops

| K index | old | new | count |
|---|---|---|---|
| 1 | R | R | 80 |
| 2 | C | C | 93 |
| 2 | R | R | 63 |
| 2 | R | U | 1 |
| 2 | U | C | 12 |
| 2 | U | R | 3 |
| 2 | U | U | 13 |
| 3 | C | C | 107 |
| 3 | R | R | 65 |
| 3 | R | U | 2 |
| 3 | U | C | 11 |
| 3 | U | R | 2 |
| 3 | U | U | 13 |
| 4 | C | C | 112 |
| 4 | R | R | 67 |
| 4 | R | U | 2 |
| 4 | U | C | 8 |
| 4 | U | R | 4 |
| 4 | U | U | 7 |
| 5 | C | C | 82 |
| 5 | R | R | 35 |
| 5 | U | C | 4 |
| 5 | U | R | 3 |
| 5 | U | U | 1 |
| 6 | C | C | 80 |

Panel: 52 orbits changed class; grid: 0. Dropped: 171 (control 13, unresolved 158); relabelling data missing: 0.

### Sets of Cal (pooled, panel, grid)

| set | pooled | panel | grid |
|---|---|---|---|
| cal | 3016 | 831 | 2185 |
| S2 | 452 | 165 | 287 |
| T_s | 31 | 16 | 15 |
| W2 | 228 | 139 | 89 |
| R_W | 197 | 123 | 74 |
| regular | 1166 | 322 | 844 |
| gate | 80 | 80 | 0 |

### Feasible set F and the selection

- elements of Theta: 2088; cost cap of F1: 185395200/1 (1.85395e8)
- elements satisfying F1 to F6 alone: 690, 2018, 2088, 2088, 1897, 1177
- |F| = 382; outcome: selected
- J_min: 0/1 (0.0); |Theta_tol| = 20
- theta* (index 371): V = SALI, FLI, WBA; qC = 1, qR = 3, w = 2, v = 0; C_R = 72430080/1 (7.24301e7)
- calibration error e = 7/1432 (0.00488827); unresolved share on Cal 19/377 (0.0503979), on S2 0/1 (0.0); gate count 0; J = 259/14136 (0.018322)


### HB analogue on Cal at theta* (C7)

Set T_s:

| X | b | cost | capped | cost-unmatched | FR_X | FR_rule | n10 | n01 | p | margin | coverage |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SALI | 51200 | 136089600/1 | false | false | 8/31 | 1/31 | 1 | 8 | 5/256 (0.0195312) | true | 25/31 |
| FLI | 51200 | 136089600/1 | false | false | 12/31 | 1/31 | 1 | 12 | 7/4096 (0.00170898) | true | 25/31 |
| MEGNO | 51200 | 136089600/1 | false | false | 7/31 | 1/31 | 1 | 7 | 9/256 (0.0351562) | true | 25/31 |
| E1 | 6400 | 107842560/1 | false | false | 25/31 | 1/31 | 0 | 24 | 1/16777216 (5.96046e-8) | true | 25/31 |
| E2 | 51200 | 92697600/1 | false | false | 23/31 | 1/31 | 0 | 22 | 1/4194304 (2.38419e-7) | true | 25/31 |
| WBA | 25600 | 76815360/1 | false | false | 4/31 | 1/31 | 1 | 4 | 3/16 (0.1875) | true | 25/31 |

Set W2:

| X | b | cost | capped | cost-unmatched | FR_X | FR_rule | n10 | n01 | p | margin | coverage |
|---|---|---|---|---|---|---|---|---|---|---|---|
| LE | 102400 | 21258240/1 | true | true | 3/38 | 1/228 | 1 | 18 | 5/131072 (3.8147e-5) | true | 109/114 |
| SALI | 51200 | 136089600/1 | false | false | 11/114 | 1/228 | 1 | 22 | 3/1048576 (2.86102e-6) | true | 109/114 |
| FLI | 51200 | 136089600/1 | false | false | 9/76 | 1/228 | 1 | 27 | 29/268435456 (1.08033e-7) | true | 109/114 |
| MEGNO | 51200 | 136089600/1 | false | false | 29/228 | 1/228 | 1 | 29 | 31/1073741824 (2.8871e-8) | true | 109/114 |
| E1 | 6400 | 107842560/1 | false | false | 23/114 | 1/228 | 0 | 45 | 1/35184372088832 (2.84217e-14) | true | 109/114 |
| E2 | 51200 | 92697600/1 | false | false | 5/38 | 1/228 | 0 | 29 | 1/536870912 (1.86265e-9) | true | 109/114 |
| WBA | 25600 | 76815360/1 | false | false | 1/38 | 1/228 | 1 | 6 | 1/16 (0.0625) | true | 109/114 |

Risk and coverage of the rule family around theta* and of the single indicators:

| qC | qR | coverage | e |
|---|---|---|---|
| 1 | 2 | 2985/3016 | 4/597 |
| 1 | 3 | 358/377 | 7/1432 |
| 2 | 3 | 95/104 | 1/2755 |

| X | b | coverage | e |
|---|---|---|---|
| LE | 102400 | 375/377 | 9/1000 |
| SALI | 51200 | 1/1 | 5/377 |
| FLI | 51200 | 1/1 | 4/377 |
| MEGNO | 51200 | 2973/3016 | 113/2973 |
| E1 | 6400 | 3011/3016 | 82/3011 |
| E2 | 51200 | 1/1 | 17/1508 |
| WBA | 25600 | 3011/3016 | 9/3011 |

First 20 elements of Theta_tol in the order of C6: 371, 407, 443, 372, 408, 444, 971, 972, 587, 659, 588, 660, 1163, 1211, 1259, 1164, 1212, 1260, 1637, 1638


### PP for HB and HC

- optimism: used 200, skipped 0; r_corr T_s 0.0591379262614, W2 0.00890288810764, R_W 0.000997882103743
- rule rates: T_s 0.0591379262614, R_W 0.01; targets: cap, cap true, n_T 285, n_W 1235

| j | rho | nT | nRW | nW | allowed | successes | runs | power | se |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1/1 | 57 | 190 | 247 | false | undefined | 20000 | undefined | undefined |
| 1 | 5/4 | 71 | 238 | 309 | false | undefined | 20000 | undefined | undefined |
| 2 | 3/2 | 86 | 285 | 371 | false | undefined | 20000 | undefined | undefined |
| 3 | 7/4 | 100 | 333 | 433 | true | 6250 | 20000 | 0.3125 | 0.0032775276505317237 |
| 4 | 2/1 | 114 | 380 | 494 | true | 6604 | 20000 | 0.3302 | 0.003325416966336703 |
| 5 | 9/4 | 128 | 428 | 556 | true | 6950 | 20000 | 0.3475 | 0.003367074026510258 |
| 6 | 5/2 | 143 | 475 | 618 | true | 7180 | 20000 | 0.359 | 0.003392042157756887 |
| 7 | 11/4 | 157 | 523 | 680 | true | 7204 | 20000 | 0.3602 | 0.003394524708998301 |
| 8 | 3/1 | 171 | 570 | 741 | true | 7037 | 20000 | 0.35185 | 0.0033767704800593125 |
| 9 | 13/4 | 185 | 618 | 803 | true | 7057 | 20000 | 0.35285 | 0.0033789560333037775 |
| 10 | 7/2 | 200 | 665 | 865 | true | 6855 | 20000 | 0.34275 | 0.0033561319811652225 |
| 11 | 15/4 | 214 | 713 | 927 | true | 6842 | 20000 | 0.3421 | 0.0033546057145363594 |
| 12 | 4/1 | 228 | 760 | 988 | true | 6629 | 20000 | 0.33145 | 0.003328595030189164 |
| 13 | 17/4 | 242 | 808 | 1050 | true | 6448 | 20000 | 0.3224 | 0.00330498290464565 |
| 14 | 9/2 | 257 | 855 | 1112 | true | 6403 | 20000 | 0.32015 | 0.0032988935834609764 |
| 15 | 19/4 | 271 | 903 | 1174 | true | 6282 | 20000 | 0.3141 | 0.0032820815803389163 |
| 16 | 5/1 | 285 | 950 | 1235 | true | 6163 | 20000 | 0.30815 | 0.0032649163656975964 |

Type-I check (joint and single rejection rates): maximum joint 0.0005

- HC: AUC_I_cal 0.9146362184388487 (failures 51, non-failures 2885), AUC_I_used 0.7073181092194243, target n_f 45, underpowered false


## F1 data

| index | theta | C_R | FR(W2; theta) | FR(T_s; theta) | J | in Theta_tol | theta* |
|---|---|---|---|---|---|---|---|
| 3 | V = SALI, FLI, MEGNO; qC = 1, qR = 2, w = 1, v = 0 | 34022400/1 | 5/76 | 5/31 | 535/4712 | false | false |
| 4 | V = SALI, FLI, MEGNO; qC = 1, qR = 2, w = 1, v = 1 | 34022400/1 | 3/38 | 7/31 | 359/2356 | false | false |
| 5 | V = SALI, FLI, MEGNO; qC = 1, qR = 2, w = 2, v = 0 | 34022400/1 | 3/76 | 4/31 | 397/4712 | false | false |
| 6 | V = SALI, FLI, MEGNO; qC = 1, qR = 2, w = 2, v = 1 | 34022400/1 | 7/114 | 8/31 | 1129/7068 | false | false |
| 9 | V = SALI, FLI, MEGNO; qC = 1, qR = 3, w = 1, v = 0 | 34022400/1 | 7/114 | 4/31 | 673/7068 | false | false |
| 10 | V = SALI, FLI, MEGNO; qC = 1, qR = 3, w = 1, v = 1 | 34022400/1 | 7/114 | 4/31 | 673/7068 | false | false |
| 11 | V = SALI, FLI, MEGNO; qC = 1, qR = 3, w = 2, v = 0 | 34022400/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 12 | V = SALI, FLI, MEGNO; qC = 1, qR = 3, w = 2, v = 1 | 34022400/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 13 | V = SALI, FLI, MEGNO; qC = 2, qR = 2, w = 0, v = 0 | 34022400/1 | 23/228 | 6/31 | 2081/14136 | false | false |
| 14 | V = SALI, FLI, MEGNO; qC = 2, qR = 2, w = 0, v = 1 | 34022400/1 | 11/76 | 13/31 | 1329/4712 | false | false |
| 15 | V = SALI, FLI, MEGNO; qC = 2, qR = 2, w = 1, v = 0 | 34022400/1 | 5/76 | 5/31 | 535/4712 | false | false |
| 16 | V = SALI, FLI, MEGNO; qC = 2, qR = 2, w = 1, v = 1 | 34022400/1 | 11/114 | 11/31 | 1595/7068 | false | false |
| 19 | V = SALI, FLI, MEGNO; qC = 2, qR = 3, w = 0, v = 0 | 34022400/1 | 23/228 | 6/31 | 2081/14136 | false | false |
| 20 | V = SALI, FLI, MEGNO; qC = 2, qR = 3, w = 0, v = 1 | 34022400/1 | 23/228 | 6/31 | 2081/14136 | false | false |
| 21 | V = SALI, FLI, MEGNO; qC = 2, qR = 3, w = 1, v = 0 | 34022400/1 | 7/114 | 4/31 | 673/7068 | false | false |
| 22 | V = SALI, FLI, MEGNO; qC = 2, qR = 3, w = 1, v = 1 | 34022400/1 | 7/114 | 4/31 | 673/7068 | false | false |
| 25 | V = SALI, FLI, MEGNO; qC = 3, qR = 2, w = 0, v = 0 | 34022400/1 | 23/228 | 6/31 | 2081/14136 | false | false |
| 26 | V = SALI, FLI, MEGNO; qC = 3, qR = 2, w = 0, v = 1 | 34022400/1 | 11/76 | 13/31 | 1329/4712 | false | false |
| 28 | V = SALI, FLI, MEGNO; qC = 3, qR = 2, w = 1, v = 1 | 34022400/1 | 11/114 | 11/31 | 1595/7068 | false | false |
| 31 | V = SALI, FLI, MEGNO; qC = 3, qR = 3, w = 0, v = 0 | 34022400/1 | 23/228 | 6/31 | 2081/14136 | false | false |
| 32 | V = SALI, FLI, MEGNO; qC = 3, qR = 3, w = 0, v = 1 | 34022400/1 | 23/228 | 6/31 | 2081/14136 | false | false |
| 145 | V = SALI, FLI, E2; qC = 1, qR = 2, w = 0, v = 0 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 146 | V = SALI, FLI, E2; qC = 1, qR = 2, w = 0, v = 1 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 147 | V = SALI, FLI, E2; qC = 1, qR = 2, w = 1, v = 0 | 57196800/1 | 17/228 | 10/31 | 2807/14136 | false | false |
| 148 | V = SALI, FLI, E2; qC = 1, qR = 2, w = 1, v = 1 | 57196800/1 | 5/57 | 13/31 | 448/1767 | false | false |
| 149 | V = SALI, FLI, E2; qC = 1, qR = 2, w = 2, v = 0 | 57196800/1 | 4/57 | 10/31 | 347/1767 | false | false |
| 150 | V = SALI, FLI, E2; qC = 1, qR = 2, w = 2, v = 1 | 57196800/1 | 1/12 | 13/31 | 187/744 | false | false |
| 151 | V = SALI, FLI, E2; qC = 1, qR = 3, w = 0, v = 0 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 152 | V = SALI, FLI, E2; qC = 1, qR = 3, w = 0, v = 1 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 153 | V = SALI, FLI, E2; qC = 1, qR = 3, w = 1, v = 0 | 57196800/1 | 7/114 | 9/31 | 1243/7068 | false | false |
| 154 | V = SALI, FLI, E2; qC = 1, qR = 3, w = 1, v = 1 | 57196800/1 | 7/114 | 9/31 | 1243/7068 | false | false |
| 155 | V = SALI, FLI, E2; qC = 1, qR = 3, w = 2, v = 0 | 57196800/1 | 1/19 | 7/31 | 82/589 | false | false |
| 156 | V = SALI, FLI, E2; qC = 1, qR = 3, w = 2, v = 1 | 57196800/1 | 1/19 | 7/31 | 82/589 | false | false |
| 157 | V = SALI, FLI, E2; qC = 2, qR = 2, w = 0, v = 0 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 158 | V = SALI, FLI, E2; qC = 2, qR = 2, w = 0, v = 1 | 57196800/1 | 3/19 | 17/31 | 208/589 | false | false |
| 159 | V = SALI, FLI, E2; qC = 2, qR = 2, w = 1, v = 0 | 57196800/1 | 17/228 | 10/31 | 2807/14136 | false | false |
| 160 | V = SALI, FLI, E2; qC = 2, qR = 2, w = 1, v = 1 | 57196800/1 | 9/76 | 15/31 | 1419/4712 | false | false |
| 161 | V = SALI, FLI, E2; qC = 2, qR = 2, w = 2, v = 0 | 57196800/1 | 4/57 | 10/31 | 347/1767 | false | false |
| 162 | V = SALI, FLI, E2; qC = 2, qR = 2, w = 2, v = 1 | 57196800/1 | 7/76 | 14/31 | 1281/4712 | false | false |
| 163 | V = SALI, FLI, E2; qC = 2, qR = 3, w = 0, v = 0 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 164 | V = SALI, FLI, E2; qC = 2, qR = 3, w = 0, v = 1 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 165 | V = SALI, FLI, E2; qC = 2, qR = 3, w = 1, v = 0 | 57196800/1 | 7/114 | 9/31 | 1243/7068 | false | false |
| 166 | V = SALI, FLI, E2; qC = 2, qR = 3, w = 1, v = 1 | 57196800/1 | 7/114 | 9/31 | 1243/7068 | false | false |
| 169 | V = SALI, FLI, E2; qC = 3, qR = 2, w = 0, v = 0 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 170 | V = SALI, FLI, E2; qC = 3, qR = 2, w = 0, v = 1 | 57196800/1 | 3/19 | 17/31 | 208/589 | false | false |
| 172 | V = SALI, FLI, E2; qC = 3, qR = 2, w = 1, v = 1 | 57196800/1 | 9/76 | 15/31 | 1419/4712 | false | false |
| 174 | V = SALI, FLI, E2; qC = 3, qR = 2, w = 2, v = 1 | 57196800/1 | 7/76 | 14/31 | 1281/4712 | false | false |
| 175 | V = SALI, FLI, E2; qC = 3, qR = 3, w = 0, v = 0 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 176 | V = SALI, FLI, E2; qC = 3, qR = 3, w = 0, v = 1 | 57196800/1 | 1/12 | 11/31 | 163/744 | false | false |
| 183 | V = SALI, MEGNO, E2; qC = 1, qR = 2, w = 1, v = 0 | 57196800/1 | 3/76 | 4/31 | 397/4712 | false | false |
| 184 | V = SALI, MEGNO, E2; qC = 1, qR = 2, w = 1, v = 1 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 185 | V = SALI, MEGNO, E2; qC = 1, qR = 2, w = 2, v = 0 | 57196800/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 186 | V = SALI, MEGNO, E2; qC = 1, qR = 2, w = 2, v = 1 | 57196800/1 | 13/228 | 8/31 | 2227/14136 | false | false |
| 189 | V = SALI, MEGNO, E2; qC = 1, qR = 3, w = 1, v = 0 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 190 | V = SALI, MEGNO, E2; qC = 1, qR = 3, w = 1, v = 1 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 191 | V = SALI, MEGNO, E2; qC = 1, qR = 3, w = 2, v = 0 | 57196800/1 | 1/57 | 1/31 | 44/1767 | false | false |
| 192 | V = SALI, MEGNO, E2; qC = 1, qR = 3, w = 2, v = 1 | 57196800/1 | 1/57 | 1/31 | 44/1767 | false | false |
| 193 | V = SALI, MEGNO, E2; qC = 2, qR = 2, w = 0, v = 0 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 194 | V = SALI, MEGNO, E2; qC = 2, qR = 2, w = 0, v = 1 | 57196800/1 | 11/76 | 13/31 | 1329/4712 | false | false |
| 195 | V = SALI, MEGNO, E2; qC = 2, qR = 2, w = 1, v = 0 | 57196800/1 | 3/76 | 4/31 | 397/4712 | false | false |
| 196 | V = SALI, MEGNO, E2; qC = 2, qR = 2, w = 1, v = 1 | 57196800/1 | 11/114 | 11/31 | 1595/7068 | false | false |
| 199 | V = SALI, MEGNO, E2; qC = 2, qR = 3, w = 0, v = 0 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 200 | V = SALI, MEGNO, E2; qC = 2, qR = 3, w = 0, v = 1 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 201 | V = SALI, MEGNO, E2; qC = 2, qR = 3, w = 1, v = 0 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 202 | V = SALI, MEGNO, E2; qC = 2, qR = 3, w = 1, v = 1 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 206 | V = SALI, MEGNO, E2; qC = 3, qR = 2, w = 0, v = 1 | 57196800/1 | 11/76 | 13/31 | 1329/4712 | false | false |
| 217 | V = FLI, MEGNO, E2; qC = 1, qR = 2, w = 0, v = 0 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 218 | V = FLI, MEGNO, E2; qC = 1, qR = 2, w = 0, v = 1 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 219 | V = FLI, MEGNO, E2; qC = 1, qR = 2, w = 1, v = 0 | 57196800/1 | 3/76 | 4/31 | 397/4712 | false | false |
| 220 | V = FLI, MEGNO, E2; qC = 1, qR = 2, w = 1, v = 1 | 57196800/1 | 5/76 | 9/31 | 839/4712 | false | false |
| 221 | V = FLI, MEGNO, E2; qC = 1, qR = 2, w = 2, v = 0 | 57196800/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 222 | V = FLI, MEGNO, E2; qC = 1, qR = 2, w = 2, v = 1 | 57196800/1 | 17/228 | 12/31 | 3263/14136 | false | false |
| 223 | V = FLI, MEGNO, E2; qC = 1, qR = 3, w = 0, v = 0 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 224 | V = FLI, MEGNO, E2; qC = 1, qR = 3, w = 0, v = 1 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 225 | V = FLI, MEGNO, E2; qC = 1, qR = 3, w = 1, v = 0 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 226 | V = FLI, MEGNO, E2; qC = 1, qR = 3, w = 1, v = 1 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 229 | V = FLI, MEGNO, E2; qC = 2, qR = 2, w = 0, v = 0 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 230 | V = FLI, MEGNO, E2; qC = 2, qR = 2, w = 0, v = 1 | 57196800/1 | 1/6 | 17/31 | 133/372 | false | false |
| 232 | V = FLI, MEGNO, E2; qC = 2, qR = 2, w = 1, v = 1 | 57196800/1 | 9/76 | 15/31 | 1419/4712 | false | false |
| 234 | V = FLI, MEGNO, E2; qC = 2, qR = 2, w = 2, v = 1 | 57196800/1 | 7/76 | 14/31 | 1281/4712 | false | false |
| 235 | V = FLI, MEGNO, E2; qC = 2, qR = 3, w = 0, v = 0 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 236 | V = FLI, MEGNO, E2; qC = 2, qR = 3, w = 0, v = 1 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 242 | V = FLI, MEGNO, E2; qC = 3, qR = 2, w = 0, v = 1 | 57196800/1 | 1/6 | 17/31 | 133/372 | false | false |
| 244 | V = FLI, MEGNO, E2; qC = 3, qR = 2, w = 1, v = 1 | 57196800/1 | 9/76 | 15/31 | 1419/4712 | false | false |
| 246 | V = FLI, MEGNO, E2; qC = 3, qR = 2, w = 2, v = 1 | 57196800/1 | 7/76 | 14/31 | 1281/4712 | false | false |
| 361 | V = SALI, FLI, WBA; qC = 1, qR = 2, w = 0, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 362 | V = SALI, FLI, WBA; qC = 1, qR = 2, w = 0, v = 1 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 363 | V = SALI, FLI, WBA; qC = 1, qR = 2, w = 1, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 364 | V = SALI, FLI, WBA; qC = 1, qR = 2, w = 1, v = 1 | 72430080/1 | 5/114 | 6/31 | 839/7068 | false | false |
| 365 | V = SALI, FLI, WBA; qC = 1, qR = 2, w = 2, v = 0 | 72430080/1 | 7/228 | 4/31 | 1129/14136 | false | false |
| 366 | V = SALI, FLI, WBA; qC = 1, qR = 2, w = 2, v = 1 | 72430080/1 | 2/57 | 5/31 | 347/3534 | false | false |
| 367 | V = SALI, FLI, WBA; qC = 1, qR = 3, w = 0, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 368 | V = SALI, FLI, WBA; qC = 1, qR = 3, w = 0, v = 1 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 369 | V = SALI, FLI, WBA; qC = 1, qR = 3, w = 1, v = 0 | 72430080/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 370 | V = SALI, FLI, WBA; qC = 1, qR = 3, w = 1, v = 1 | 72430080/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 371 | V = SALI, FLI, WBA; qC = 1, qR = 3, w = 2, v = 0 | 72430080/1 | 1/228 | 1/31 | 259/14136 | true | true |
| 372 | V = SALI, FLI, WBA; qC = 1, qR = 3, w = 2, v = 1 | 72430080/1 | 1/228 | 1/31 | 259/14136 | true | false |
| 373 | V = SALI, FLI, WBA; qC = 2, qR = 2, w = 0, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 374 | V = SALI, FLI, WBA; qC = 2, qR = 2, w = 0, v = 1 | 72430080/1 | 17/114 | 15/31 | 2237/7068 | false | false |
| 375 | V = SALI, FLI, WBA; qC = 2, qR = 2, w = 1, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 376 | V = SALI, FLI, WBA; qC = 2, qR = 2, w = 1, v = 1 | 72430080/1 | 25/228 | 13/31 | 3739/14136 | false | false |
| 378 | V = SALI, FLI, WBA; qC = 2, qR = 2, w = 2, v = 1 | 72430080/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 379 | V = SALI, FLI, WBA; qC = 2, qR = 3, w = 0, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 380 | V = SALI, FLI, WBA; qC = 2, qR = 3, w = 0, v = 1 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 381 | V = SALI, FLI, WBA; qC = 2, qR = 3, w = 1, v = 0 | 72430080/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 382 | V = SALI, FLI, WBA; qC = 2, qR = 3, w = 1, v = 1 | 72430080/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 385 | V = SALI, FLI, WBA; qC = 3, qR = 2, w = 0, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 386 | V = SALI, FLI, WBA; qC = 3, qR = 2, w = 0, v = 1 | 72430080/1 | 17/114 | 15/31 | 2237/7068 | false | false |
| 388 | V = SALI, FLI, WBA; qC = 3, qR = 2, w = 1, v = 1 | 72430080/1 | 25/228 | 13/31 | 3739/14136 | false | false |
| 390 | V = SALI, FLI, WBA; qC = 3, qR = 2, w = 2, v = 1 | 72430080/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 391 | V = SALI, FLI, WBA; qC = 3, qR = 3, w = 0, v = 0 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 392 | V = SALI, FLI, WBA; qC = 3, qR = 3, w = 0, v = 1 | 72430080/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 399 | V = SALI, MEGNO, WBA; qC = 1, qR = 2, w = 1, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 400 | V = SALI, MEGNO, WBA; qC = 1, qR = 2, w = 1, v = 1 | 72430080/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 401 | V = SALI, MEGNO, WBA; qC = 1, qR = 2, w = 2, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 402 | V = SALI, MEGNO, WBA; qC = 1, qR = 2, w = 2, v = 1 | 72430080/1 | 1/57 | 2/31 | 145/3534 | false | false |
| 405 | V = SALI, MEGNO, WBA; qC = 1, qR = 3, w = 1, v = 0 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 406 | V = SALI, MEGNO, WBA; qC = 1, qR = 3, w = 1, v = 1 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 407 | V = SALI, MEGNO, WBA; qC = 1, qR = 3, w = 2, v = 0 | 72430080/1 | 0/1 | 0/1 | 0/1 | true | false |
| 408 | V = SALI, MEGNO, WBA; qC = 1, qR = 3, w = 2, v = 1 | 72430080/1 | 0/1 | 0/1 | 0/1 | true | false |
| 409 | V = SALI, MEGNO, WBA; qC = 2, qR = 2, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 410 | V = SALI, MEGNO, WBA; qC = 2, qR = 2, w = 0, v = 1 | 72430080/1 | 7/57 | 9/31 | 365/1767 | false | false |
| 411 | V = SALI, MEGNO, WBA; qC = 2, qR = 2, w = 1, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 412 | V = SALI, MEGNO, WBA; qC = 2, qR = 2, w = 1, v = 1 | 72430080/1 | 3/38 | 7/31 | 359/2356 | false | false |
| 414 | V = SALI, MEGNO, WBA; qC = 2, qR = 2, w = 2, v = 1 | 72430080/1 | 7/228 | 3/31 | 901/14136 | false | false |
| 415 | V = SALI, MEGNO, WBA; qC = 2, qR = 3, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 416 | V = SALI, MEGNO, WBA; qC = 2, qR = 3, w = 0, v = 1 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 417 | V = SALI, MEGNO, WBA; qC = 2, qR = 3, w = 1, v = 0 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 418 | V = SALI, MEGNO, WBA; qC = 2, qR = 3, w = 1, v = 1 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 422 | V = SALI, MEGNO, WBA; qC = 3, qR = 2, w = 0, v = 1 | 72430080/1 | 7/57 | 9/31 | 365/1767 | false | false |
| 424 | V = SALI, MEGNO, WBA; qC = 3, qR = 2, w = 1, v = 1 | 72430080/1 | 3/38 | 7/31 | 359/2356 | false | false |
| 433 | V = FLI, MEGNO, WBA; qC = 1, qR = 2, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 434 | V = FLI, MEGNO, WBA; qC = 1, qR = 2, w = 0, v = 1 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 435 | V = FLI, MEGNO, WBA; qC = 1, qR = 2, w = 1, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 436 | V = FLI, MEGNO, WBA; qC = 1, qR = 2, w = 1, v = 1 | 72430080/1 | 7/228 | 4/31 | 1129/14136 | false | false |
| 437 | V = FLI, MEGNO, WBA; qC = 1, qR = 2, w = 2, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 438 | V = FLI, MEGNO, WBA; qC = 1, qR = 2, w = 2, v = 1 | 72430080/1 | 5/228 | 3/31 | 839/14136 | false | false |
| 439 | V = FLI, MEGNO, WBA; qC = 1, qR = 3, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 440 | V = FLI, MEGNO, WBA; qC = 1, qR = 3, w = 0, v = 1 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 441 | V = FLI, MEGNO, WBA; qC = 1, qR = 3, w = 1, v = 0 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 442 | V = FLI, MEGNO, WBA; qC = 1, qR = 3, w = 1, v = 1 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 443 | V = FLI, MEGNO, WBA; qC = 1, qR = 3, w = 2, v = 0 | 72430080/1 | 0/1 | 0/1 | 0/1 | true | false |
| 444 | V = FLI, MEGNO, WBA; qC = 1, qR = 3, w = 2, v = 1 | 72430080/1 | 0/1 | 0/1 | 0/1 | true | false |
| 445 | V = FLI, MEGNO, WBA; qC = 2, qR = 2, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 446 | V = FLI, MEGNO, WBA; qC = 2, qR = 2, w = 0, v = 1 | 72430080/1 | 11/76 | 13/31 | 1329/4712 | false | false |
| 447 | V = FLI, MEGNO, WBA; qC = 2, qR = 2, w = 1, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 448 | V = FLI, MEGNO, WBA; qC = 2, qR = 2, w = 1, v = 1 | 72430080/1 | 7/76 | 9/31 | 901/4712 | false | false |
| 450 | V = FLI, MEGNO, WBA; qC = 2, qR = 2, w = 2, v = 1 | 72430080/1 | 5/114 | 6/31 | 839/7068 | false | false |
| 451 | V = FLI, MEGNO, WBA; qC = 2, qR = 3, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 452 | V = FLI, MEGNO, WBA; qC = 2, qR = 3, w = 0, v = 1 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 453 | V = FLI, MEGNO, WBA; qC = 2, qR = 3, w = 1, v = 0 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 454 | V = FLI, MEGNO, WBA; qC = 2, qR = 3, w = 1, v = 1 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 458 | V = FLI, MEGNO, WBA; qC = 3, qR = 2, w = 0, v = 1 | 72430080/1 | 11/76 | 13/31 | 1329/4712 | false | false |
| 460 | V = FLI, MEGNO, WBA; qC = 3, qR = 2, w = 1, v = 1 | 72430080/1 | 7/76 | 9/31 | 901/4712 | false | false |
| 577 | V = SALI, E2, WBA; qC = 1, qR = 2, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 578 | V = SALI, E2, WBA; qC = 1, qR = 2, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 579 | V = SALI, E2, WBA; qC = 1, qR = 2, w = 1, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 580 | V = SALI, E2, WBA; qC = 1, qR = 2, w = 1, v = 1 | 95604480/1 | 5/114 | 6/31 | 839/7068 | false | false |
| 581 | V = SALI, E2, WBA; qC = 1, qR = 2, w = 2, v = 0 | 95604480/1 | 7/228 | 4/31 | 1129/14136 | false | false |
| 582 | V = SALI, E2, WBA; qC = 1, qR = 2, w = 2, v = 1 | 95604480/1 | 2/57 | 5/31 | 347/3534 | false | false |
| 583 | V = SALI, E2, WBA; qC = 1, qR = 3, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 584 | V = SALI, E2, WBA; qC = 1, qR = 3, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 585 | V = SALI, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 586 | V = SALI, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 587 | V = SALI, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 1/228 | 1/31 | 259/14136 | true | false |
| 588 | V = SALI, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 1/228 | 1/31 | 259/14136 | true | false |
| 589 | V = SALI, E2, WBA; qC = 2, qR = 2, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 590 | V = SALI, E2, WBA; qC = 2, qR = 2, w = 0, v = 1 | 95604480/1 | 11/114 | 14/31 | 1937/7068 | false | false |
| 591 | V = SALI, E2, WBA; qC = 2, qR = 2, w = 1, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 592 | V = SALI, E2, WBA; qC = 2, qR = 2, w = 1, v = 1 | 95604480/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 594 | V = SALI, E2, WBA; qC = 2, qR = 2, w = 2, v = 1 | 95604480/1 | 5/76 | 10/31 | 915/4712 | false | false |
| 595 | V = SALI, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 596 | V = SALI, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 597 | V = SALI, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 598 | V = SALI, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 602 | V = SALI, E2, WBA; qC = 3, qR = 2, w = 0, v = 1 | 95604480/1 | 11/114 | 14/31 | 1937/7068 | false | false |
| 604 | V = SALI, E2, WBA; qC = 3, qR = 2, w = 1, v = 1 | 95604480/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 613 | V = FLI, E2, WBA; qC = 1, qR = 2, w = 0, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 614 | V = FLI, E2, WBA; qC = 1, qR = 2, w = 0, v = 1 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 615 | V = FLI, E2, WBA; qC = 1, qR = 2, w = 1, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 616 | V = FLI, E2, WBA; qC = 1, qR = 2, w = 1, v = 1 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 617 | V = FLI, E2, WBA; qC = 1, qR = 2, w = 2, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 618 | V = FLI, E2, WBA; qC = 1, qR = 2, w = 2, v = 1 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 619 | V = FLI, E2, WBA; qC = 1, qR = 3, w = 0, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 620 | V = FLI, E2, WBA; qC = 1, qR = 3, w = 0, v = 1 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 621 | V = FLI, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 3/76 | 6/31 | 549/4712 | false | false |
| 622 | V = FLI, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 3/76 | 6/31 | 549/4712 | false | false |
| 623 | V = FLI, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 1/57 | 4/31 | 259/3534 | false | false |
| 624 | V = FLI, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 1/57 | 4/31 | 259/3534 | false | false |
| 625 | V = FLI, E2, WBA; qC = 2, qR = 2, w = 0, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 626 | V = FLI, E2, WBA; qC = 2, qR = 2, w = 0, v = 1 | 95604480/1 | 2/19 | 16/31 | 183/589 | false | false |
| 627 | V = FLI, E2, WBA; qC = 2, qR = 2, w = 1, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 628 | V = FLI, E2, WBA; qC = 2, qR = 2, w = 1, v = 1 | 95604480/1 | 1/12 | 13/31 | 187/744 | false | false |
| 630 | V = FLI, E2, WBA; qC = 2, qR = 2, w = 2, v = 1 | 95604480/1 | 1/12 | 13/31 | 187/744 | false | false |
| 631 | V = FLI, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 632 | V = FLI, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 11/228 | 7/31 | 1937/14136 | false | false |
| 633 | V = FLI, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 3/76 | 6/31 | 549/4712 | false | false |
| 634 | V = FLI, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 3/76 | 6/31 | 549/4712 | false | false |
| 638 | V = FLI, E2, WBA; qC = 3, qR = 2, w = 0, v = 1 | 95604480/1 | 2/19 | 16/31 | 183/589 | false | false |
| 640 | V = FLI, E2, WBA; qC = 3, qR = 2, w = 1, v = 1 | 95604480/1 | 1/12 | 13/31 | 187/744 | false | false |
| 642 | V = FLI, E2, WBA; qC = 3, qR = 2, w = 2, v = 1 | 95604480/1 | 1/12 | 13/31 | 187/744 | false | false |
| 649 | V = MEGNO, E2, WBA; qC = 1, qR = 2, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 650 | V = MEGNO, E2, WBA; qC = 1, qR = 2, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 651 | V = MEGNO, E2, WBA; qC = 1, qR = 2, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 652 | V = MEGNO, E2, WBA; qC = 1, qR = 2, w = 1, v = 1 | 95604480/1 | 7/228 | 4/31 | 1129/14136 | false | false |
| 653 | V = MEGNO, E2, WBA; qC = 1, qR = 2, w = 2, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 654 | V = MEGNO, E2, WBA; qC = 1, qR = 2, w = 2, v = 1 | 95604480/1 | 5/228 | 3/31 | 839/14136 | false | false |
| 655 | V = MEGNO, E2, WBA; qC = 1, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 656 | V = MEGNO, E2, WBA; qC = 1, qR = 3, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 657 | V = MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 658 | V = MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 659 | V = MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 660 | V = MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 661 | V = MEGNO, E2, WBA; qC = 2, qR = 2, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 662 | V = MEGNO, E2, WBA; qC = 2, qR = 2, w = 0, v = 1 | 95604480/1 | 7/76 | 12/31 | 1129/4712 | false | false |
| 663 | V = MEGNO, E2, WBA; qC = 2, qR = 2, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 664 | V = MEGNO, E2, WBA; qC = 2, qR = 2, w = 1, v = 1 | 95604480/1 | 13/228 | 7/31 | 1999/14136 | false | false |
| 666 | V = MEGNO, E2, WBA; qC = 2, qR = 2, w = 2, v = 1 | 95604480/1 | 2/57 | 5/31 | 347/3534 | false | false |
| 667 | V = MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 668 | V = MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 669 | V = MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 670 | V = MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 674 | V = MEGNO, E2, WBA; qC = 3, qR = 2, w = 0, v = 1 | 95604480/1 | 7/76 | 12/31 | 1129/4712 | false | false |
| 771 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 3, w = 1, v = 0 | 57196800/1 | 3/76 | 4/31 | 397/4712 | false | false |
| 772 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 3, w = 1, v = 1 | 57196800/1 | 1/19 | 6/31 | 145/1178 | false | false |
| 773 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 3, w = 2, v = 0 | 57196800/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 774 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 3, w = 2, v = 1 | 57196800/1 | 13/228 | 8/31 | 2227/14136 | false | false |
| 777 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 4, w = 1, v = 0 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 778 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 4, w = 1, v = 1 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 779 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 4, w = 2, v = 0 | 57196800/1 | 1/57 | 1/31 | 44/1767 | false | false |
| 780 | V = SALI, FLI, MEGNO, E2; qC = 1, qR = 4, w = 2, v = 1 | 57196800/1 | 1/57 | 1/31 | 44/1767 | false | false |
| 781 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 3, w = 0, v = 0 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 782 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 3, w = 0, v = 1 | 57196800/1 | 8/57 | 13/31 | 989/3534 | false | false |
| 783 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 3, w = 1, v = 0 | 57196800/1 | 3/76 | 4/31 | 397/4712 | false | false |
| 784 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 3, w = 1, v = 1 | 57196800/1 | 11/114 | 11/31 | 1595/7068 | false | false |
| 787 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 4, w = 0, v = 0 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 788 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 4, w = 0, v = 1 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 789 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 4, w = 1, v = 0 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 790 | V = SALI, FLI, MEGNO, E2; qC = 2, qR = 4, w = 1, v = 1 | 57196800/1 | 1/38 | 2/31 | 107/2356 | false | false |
| 793 | V = SALI, FLI, MEGNO, E2; qC = 3, qR = 3, w = 0, v = 0 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 794 | V = SALI, FLI, MEGNO, E2; qC = 3, qR = 3, w = 0, v = 1 | 57196800/1 | 8/57 | 13/31 | 989/3534 | false | false |
| 796 | V = SALI, FLI, MEGNO, E2; qC = 3, qR = 3, w = 1, v = 1 | 57196800/1 | 11/114 | 11/31 | 1595/7068 | false | false |
| 799 | V = SALI, FLI, MEGNO, E2; qC = 3, qR = 4, w = 0, v = 0 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 800 | V = SALI, FLI, MEGNO, E2; qC = 3, qR = 4, w = 0, v = 1 | 57196800/1 | 11/228 | 5/31 | 1481/14136 | false | false |
| 806 | V = SALI, FLI, MEGNO, E2; qC = 4, qR = 3, w = 0, v = 1 | 57196800/1 | 8/57 | 13/31 | 989/3534 | false | false |
| 963 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 3, w = 1, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 964 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 3, w = 1, v = 1 | 72430080/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 965 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 3, w = 2, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 966 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 3, w = 2, v = 1 | 72430080/1 | 1/57 | 2/31 | 145/3534 | false | false |
| 969 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 4, w = 1, v = 0 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 970 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 4, w = 1, v = 1 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 971 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 4, w = 2, v = 0 | 72430080/1 | 0/1 | 0/1 | 0/1 | true | false |
| 972 | V = SALI, FLI, MEGNO, WBA; qC = 1, qR = 4, w = 2, v = 1 | 72430080/1 | 0/1 | 0/1 | 0/1 | true | false |
| 973 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 3, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 974 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 3, w = 0, v = 1 | 72430080/1 | 7/57 | 9/31 | 365/1767 | false | false |
| 975 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 3, w = 1, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 976 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 3, w = 1, v = 1 | 72430080/1 | 3/38 | 7/31 | 359/2356 | false | false |
| 978 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 3, w = 2, v = 1 | 72430080/1 | 7/228 | 3/31 | 901/14136 | false | false |
| 979 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 4, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 980 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 4, w = 0, v = 1 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 981 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 4, w = 1, v = 0 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 982 | V = SALI, FLI, MEGNO, WBA; qC = 2, qR = 4, w = 1, v = 1 | 72430080/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 985 | V = SALI, FLI, MEGNO, WBA; qC = 3, qR = 3, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 986 | V = SALI, FLI, MEGNO, WBA; qC = 3, qR = 3, w = 0, v = 1 | 72430080/1 | 7/57 | 9/31 | 365/1767 | false | false |
| 988 | V = SALI, FLI, MEGNO, WBA; qC = 3, qR = 3, w = 1, v = 1 | 72430080/1 | 3/38 | 7/31 | 359/2356 | false | false |
| 991 | V = SALI, FLI, MEGNO, WBA; qC = 3, qR = 4, w = 0, v = 0 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 992 | V = SALI, FLI, MEGNO, WBA; qC = 3, qR = 4, w = 0, v = 1 | 72430080/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 998 | V = SALI, FLI, MEGNO, WBA; qC = 4, qR = 3, w = 0, v = 1 | 72430080/1 | 7/57 | 9/31 | 365/1767 | false | false |
| 1000 | V = SALI, FLI, MEGNO, WBA; qC = 4, qR = 3, w = 1, v = 1 | 72430080/1 | 3/38 | 7/31 | 359/2356 | false | false |
| 1153 | V = SALI, FLI, E2, WBA; qC = 1, qR = 3, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1154 | V = SALI, FLI, E2, WBA; qC = 1, qR = 3, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1155 | V = SALI, FLI, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1156 | V = SALI, FLI, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 5/114 | 6/31 | 839/7068 | false | false |
| 1157 | V = SALI, FLI, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 7/228 | 4/31 | 1129/14136 | false | false |
| 1158 | V = SALI, FLI, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 2/57 | 5/31 | 347/3534 | false | false |
| 1159 | V = SALI, FLI, E2, WBA; qC = 1, qR = 4, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1160 | V = SALI, FLI, E2, WBA; qC = 1, qR = 4, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1161 | V = SALI, FLI, E2, WBA; qC = 1, qR = 4, w = 1, v = 0 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 1162 | V = SALI, FLI, E2, WBA; qC = 1, qR = 4, w = 1, v = 1 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 1163 | V = SALI, FLI, E2, WBA; qC = 1, qR = 4, w = 2, v = 0 | 95604480/1 | 1/228 | 1/31 | 259/14136 | true | false |
| 1164 | V = SALI, FLI, E2, WBA; qC = 1, qR = 4, w = 2, v = 1 | 95604480/1 | 1/228 | 1/31 | 259/14136 | true | false |
| 1165 | V = SALI, FLI, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1166 | V = SALI, FLI, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 11/114 | 14/31 | 1937/7068 | false | false |
| 1167 | V = SALI, FLI, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1168 | V = SALI, FLI, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 1170 | V = SALI, FLI, E2, WBA; qC = 2, qR = 3, w = 2, v = 1 | 95604480/1 | 5/76 | 10/31 | 915/4712 | false | false |
| 1171 | V = SALI, FLI, E2, WBA; qC = 2, qR = 4, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1172 | V = SALI, FLI, E2, WBA; qC = 2, qR = 4, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1173 | V = SALI, FLI, E2, WBA; qC = 2, qR = 4, w = 1, v = 0 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 1174 | V = SALI, FLI, E2, WBA; qC = 2, qR = 4, w = 1, v = 1 | 95604480/1 | 1/38 | 4/31 | 183/2356 | false | false |
| 1177 | V = SALI, FLI, E2, WBA; qC = 3, qR = 3, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1178 | V = SALI, FLI, E2, WBA; qC = 3, qR = 3, w = 0, v = 1 | 95604480/1 | 11/114 | 14/31 | 1937/7068 | false | false |
| 1180 | V = SALI, FLI, E2, WBA; qC = 3, qR = 3, w = 1, v = 1 | 95604480/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 1182 | V = SALI, FLI, E2, WBA; qC = 3, qR = 3, w = 2, v = 1 | 95604480/1 | 5/76 | 10/31 | 915/4712 | false | false |
| 1183 | V = SALI, FLI, E2, WBA; qC = 3, qR = 4, w = 0, v = 0 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1184 | V = SALI, FLI, E2, WBA; qC = 3, qR = 4, w = 0, v = 1 | 95604480/1 | 2/57 | 4/31 | 145/1767 | false | false |
| 1190 | V = SALI, FLI, E2, WBA; qC = 4, qR = 3, w = 0, v = 1 | 95604480/1 | 11/114 | 14/31 | 1937/7068 | false | false |
| 1192 | V = SALI, FLI, E2, WBA; qC = 4, qR = 3, w = 1, v = 1 | 95604480/1 | 17/228 | 11/31 | 3035/14136 | false | false |
| 1203 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1204 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 1205 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1206 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 1/57 | 2/31 | 145/3534 | false | false |
| 1209 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1210 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1211 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 2, v = 0 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 1212 | V = SALI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 2, v = 1 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 1213 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1214 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1215 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1216 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 5/114 | 5/31 | 725/7068 | false | false |
| 1218 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 2, v = 1 | 95604480/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 1219 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1220 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1221 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1222 | V = SALI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1225 | V = SALI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1226 | V = SALI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1228 | V = SALI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 1, v = 1 | 95604480/1 | 5/114 | 5/31 | 725/7068 | false | false |
| 1231 | V = SALI, MEGNO, E2, WBA; qC = 3, qR = 4, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1232 | V = SALI, MEGNO, E2, WBA; qC = 3, qR = 4, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1249 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1250 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1251 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1252 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 7/228 | 4/31 | 1129/14136 | false | false |
| 1253 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1254 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 5/228 | 3/31 | 839/14136 | false | false |
| 1255 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1256 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1257 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1258 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1259 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 2, v = 0 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 1260 | V = FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 2, v = 1 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 1261 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1262 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 5/57 | 12/31 | 839/3534 | false | false |
| 1263 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1264 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 13/228 | 7/31 | 1999/14136 | false | false |
| 1266 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 2, v = 1 | 95604480/1 | 2/57 | 5/31 | 347/3534 | false | false |
| 1267 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1268 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1269 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1270 | V = FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1274 | V = FLI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 0, v = 1 | 95604480/1 | 5/57 | 12/31 | 839/3534 | false | false |
| 1276 | V = FLI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 1, v = 1 | 95604480/1 | 13/228 | 7/31 | 1999/14136 | false | false |
| 1286 | V = FLI, MEGNO, E2, WBA; qC = 4, qR = 3, w = 0, v = 1 | 95604480/1 | 5/57 | 12/31 | 839/3534 | false | false |
| 1623 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1624 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 1, v = 1 | 95604480/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 1625 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1626 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 3, w = 2, v = 1 | 95604480/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 1629 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1630 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 1, v = 1 | 95604480/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 1631 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 2, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1632 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 4, w = 2, v = 1 | 95604480/1 | 1/57 | 2/31 | 145/3534 | false | false |
| 1635 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 5, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1636 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 5, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1637 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 5, w = 2, v = 0 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 1638 | V = SALI, FLI, MEGNO, E2, WBA; qC = 1, qR = 5, w = 2, v = 1 | 95604480/1 | 0/1 | 0/1 | 0/1 | true | false |
| 1639 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1640 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1641 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1642 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 1, v = 1 | 95604480/1 | 7/114 | 7/31 | 1015/7068 | false | false |
| 1644 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 3, w = 2, v = 1 | 95604480/1 | 1/19 | 7/31 | 82/589 | false | false |
| 1645 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1646 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1647 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 1, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1648 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 1, v = 1 | 95604480/1 | 5/114 | 5/31 | 725/7068 | false | false |
| 1650 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 4, w = 2, v = 1 | 95604480/1 | 5/228 | 2/31 | 611/14136 | false | false |
| 1651 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 5, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1652 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 5, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1653 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 5, w = 1, v = 0 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1654 | V = SALI, FLI, MEGNO, E2, WBA; qC = 2, qR = 5, w = 1, v = 1 | 95604480/1 | 1/114 | 1/31 | 145/7068 | false | false |
| 1657 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1658 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1660 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 3, w = 1, v = 1 | 95604480/1 | 7/114 | 7/31 | 1015/7068 | false | false |
| 1663 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 4, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1664 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 4, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1666 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 4, w = 1, v = 1 | 95604480/1 | 5/114 | 5/31 | 725/7068 | false | false |
| 1669 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 5, w = 0, v = 0 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1670 | V = SALI, FLI, MEGNO, E2, WBA; qC = 3, qR = 5, w = 0, v = 1 | 95604480/1 | 1/76 | 1/31 | 107/4712 | false | false |
| 1676 | V = SALI, FLI, MEGNO, E2, WBA; qC = 4, qR = 3, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1678 | V = SALI, FLI, MEGNO, E2, WBA; qC = 4, qR = 3, w = 1, v = 1 | 95604480/1 | 7/114 | 7/31 | 1015/7068 | false | false |
| 1682 | V = SALI, FLI, MEGNO, E2, WBA; qC = 4, qR = 4, w = 0, v = 1 | 95604480/1 | 4/57 | 8/31 | 290/1767 | false | false |
| 1684 | V = SALI, FLI, MEGNO, E2, WBA; qC = 4, qR = 4, w = 1, v = 1 | 95604480/1 | 5/114 | 5/31 | 725/7068 | false | false |

digest: 4c0ffca7bcb0bb08830831be821e661b5de6ca0a53793d45c025e773b756e0fa
