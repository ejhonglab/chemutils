# Vapor pressure cache changes (2026-10-23, estimates)

Cached `inchi2vapor_pressure` values (`inchi` -> `vapor_pressure` entries in
`~/.chemutils_cache.sqlite`, on the machine this was run on) were re-computed with the
lookup from commit fda8e3b, after those described in
[2026-10-23_vapor_pressure.md](2026-10-23_vapor_pressure.md). Compared to the lookup used
then:
- PubChem values marked as estimates (e.g. '/est/', 'estimated', 'calculated', other
  than calculations from fits to measurements) are only used if there are no measured
  values. The new 'vapor_pressure_measured' output says whether the value is from
  measurements.
- Extrapolated values are not used regardless of case (e.g. '/Extrapolated/' was used
  before).

Any vapor pressures (or anything computed from them) from chemutils before this may be
affected for the compounds listed below.

## When

2026-10-23 14:29-14:33 US Pacific (PDT, UTC-7), with
`scripts/fix_vapor_pressure_cache.py --check values --write`, on all 1049 non-null,
non-manual cached values. Backup of the cache before this:
`~/.chemutils_cache.sqlite.20261023_142933.bak`

- 31 values changed and were updated.
- 4 were deleted, because there were no longer any values to use.
- 1014 only had 'vapor_pressure_measured' added (same value).

Then at 14:34, with `--check manual --write`, the 6 remaining manually entered values
(geraniol, nerol, 2,3-butanedione, hexyl acetate, phenylacetaldehyde, and ethyl
benzoate; CIDs 637566, 643820, 650, 8908, 998, and 7165), which were the same as lookup
values, also only had 'vapor_pressure_measured' added, and are now saved as non-manual.
Backup of the cache before this: `~/.chemutils_cache.sqlite.20261023_143428.bak`

Some of the changes may also be from PubChem values that changed since the earlier
re-computation (e.g. ethanol's small change).

Of the 31 changed values, the median change was 0.08 log10 units, 10% changed by more
than 2x, and 10% by more than 10x.

## Changed entries

`name` is the name(s) previously looked up (e.g. odor names) that resolved to this
compound, or PubChem's name if none were. Vapor pressures are in kPa. `measured` is
whether the new value is from measurements (rather than only estimates).

| name | CID | old | new | new / old | measured | action | time |
|---|---|---|---|---|---|---|---|
| 1,2-ethanediol | 174 | 0.0123 | 0.0065 | 0.53 | yes | update | 2026-10-23 14:30:39 |
| 1,2-ethanedithiol | 10902 | 0.753 | 0.757 | 1.01 | yes | update | 2026-10-23 14:30:40 |
| 1,2-propanediol | 1030 | 0.0173 | 0.0106 | 0.612 | yes | update | 2026-10-23 14:30:53 |
| 1-octanethiol | 8144 | 0.0577 | 0.0583 | 1.01 | yes | update | 2026-10-23 14:32:53 |
| 1-propanethiol | 7848 | 20.6 | 20.6 | 1 | yes | update | 2026-10-23 14:30:54 |
| 1-propene | 8252 | 1.16e+03 | 1.16e+03 | 1 | yes | update | 2026-10-23 14:30:45 |
| 1-tridecene | 17095 | 0.00851 | 0.008 | 0.94 | yes | update | 2026-10-23 14:30:13 |
| 2-mercaptoethanol | 1567 | 0.234 | 0.132 | 0.562 | yes | update | 2026-10-23 14:30:39 |
| 2-methyl-2-propanethiol | 6387 | 24.1 | 19 | 0.787 | yes | update | 2026-10-23 14:30:59 |
| 2-naphthol | 8663 | 0.00102 | 0.002 | 1.96 | yes | update | 2026-10-23 14:29:54 |
| 4-biphenylamine | 7102 | 5.98e-05 | 7.7e-05 | 1.29 | yes | update | 2026-10-23 14:30:03 |
| 5-allyl-1,3-benzodioxole | 5144 | 0.00941 | 0.008 | 0.85 | yes | update | 2026-10-23 14:29:34 |
| acetic acid | 176 | 2.09 | 1.51 | 0.721 | yes | update | 2026-10-23 14:30:35 |
| allyl propyl disulfide | 16591 | 0.052 | 0.05 | 0.962 | yes | update | 2026-10-23 14:31:50 |
| anthracene | 8418 | 4.04e-05 | 8e-05 | 1.98 | yes | update | 2026-10-23 14:30:13 |
| benzo[ghi]perylene | 9117 | 0.00667 | 1.33e-10 | 2e-08 | yes | update | 2026-10-23 14:30:28 |
| caffeine | 2519 | 1.2e-07 | 1.33e-09 | 0.0111 | yes | update | 2026-10-23 14:32:34 |
| cyclopentadiene | 7612 | 45.1 | 54.5 | 1.21 | yes | update | 2026-10-23 14:31:31 |
| dibenzo[a,h]anthracene | 5889 | 0.127 | 1.33e-11 | 1.05e-10 | no (estimates only) | update | 2026-10-23 14:30:29 |
| ethanethiol | 6343 | 70.5 | 58.9 | 0.835 | yes | update | 2026-10-23 14:30:40 |
| ethanol | 702 | 7.91 | 7.9 | 0.999 | yes | update | 2026-10-23 14:30:38 |
| ethene | 6325 | 6.89e+03 | 6.83e+03 | 0.991 | yes | update | 2026-10-23 14:30:33 |
| folic acid | 135398658 | 0.00827 |  |  |  | delete | 2026-10-23 14:30:25 |
| formic acid; methanoic acid | 284 | 5.68 | 4.6 | 0.81 | yes | update | 2026-10-23 14:33:20 |
| hentriacontane | 12410 | 0.0187 |  |  |  | delete | 2026-10-23 14:30:42 |
| heptane | 8900 | 6.13 | 5.07 | 0.826 | yes | update | 2026-10-23 14:32:17 |
| hexadienal; 2,4-hexadienal | 637564 | 0.64 | 0.641 | 1 | yes | update | 2026-10-23 14:32:07 |
| l-glutamic acid | 33032 | 2.27e-09 |  |  |  | delete | 2026-10-23 14:31:39 |
| methyl 3-oxobutanoate | 7757 | 0.119 | 0.17 | 1.43 | yes | update | 2026-10-23 14:31:37 |
| methyl anthranilate | 8635 | 0.00361 | 0.00267 | 0.738 | yes | update | 2026-10-23 14:33:02 |
| methyl formate; formic acid, methyl ester | 7865 | 78.1 | 64 | 0.82 | yes | update | 2026-10-23 14:30:36 |
| n-butyl-n-nitroso-1-butanamine | 13542 | 0.00625 | 0.00533 | 0.853 | yes | update | 2026-10-23 14:32:51 |
| nonane | 8141 | 0.592 | 0.59 | 0.997 | yes | update | 2026-10-23 14:33:13 |
| phenylacetonitrile | 8794 | 0.0119 | 0.0107 | 0.899 | yes | update | 2026-10-23 14:32:56 |
| tritriacontane | 12411 | 0.0536 |  |  |  | delete | 2026-10-23 14:30:42 |
