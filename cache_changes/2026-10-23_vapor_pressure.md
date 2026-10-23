# Vapor pressure cache changes (2026-10-23)

Cached `inchi2vapor_pressure` values (`inchi` -> `vapor_pressure` entries in
`~/.chemutils_cache.sqlite`, on the machine this was run on) were re-computed with the
lookup from commit 95a66c6, which only uses PubChem values at or near 25 C. Before that,
the lookup averaged all PubChem values for a compound, measured at any temperature, so a
single value measured at a high temperature could make the result orders of magnitude
too high. See the `inchi2vapor_pressure` and `_pubchemprops_temp_rank` docstrings for
which values are used now.

Any vapor pressures (or anything computed from them, e.g. headspace concentrations) from
chemutils before this date may be affected for the compounds listed below.

## When

All times are US Pacific (PDT, UTC-7), on 2026-10-23, with
`scripts/fix_vapor_pressure_cache.py`:

1. 12:16-12:27, `--write` (before `--check` was added): all 1066 non-null, non-manual
   values were re-computed. 503 changed and were updated, and
   72 were deleted, because PubChem no longer has any vapor
   pressure values that the new lookup would use for them (for about 40 of these,
   PubChem no longer has any vapor pressure data at all, so it is not known which
   temperatures their old values came from). Backup of the cache before this:
   `~/.chemutils_cache.sqlite.20261023_121653.bak`
2. 12:34-13:06, `--check manual nulls --write`: manually entered values (from
   `add_properties_to_cache`, which had all been copied from earlier PubChem lookups)
   were re-computed too. 48 changed and were updated (now saved as non-manual), and
   1 deleted. All 6558 cached null values were also looked up again, and 7
   now have values (from values PubChem has added, or that the old lookup failed to
   parse). Backup of the cache before this: `~/.chemutils_cache.sqlite.20261023_123404.bak`
3. 13:07, `--check nulls --inchis <failed lookups> --write`: the 37 null values whose lookups
   failed in run 2 (PubChem rate limiting) were looked up again. All are still null.
   Backup of the cache before this: `~/.chemutils_cache.sqlite.20261023_130722.bak`

Of the 551 values that changed (and were not deleted), the median change was
0.07 log10 units, 21% changed by more than 2x, and
13% by more than 10x (49% decreased).

## Changed entries

`name` is the name(s) previously looked up (e.g. odor names) that resolved to this
compound, or PubChem's name if none were. `entry` is what the entry was before the
change (a non-manual value, a manually entered value, or null). Vapor pressures are in
kPa. `time` is when the entry was changed (only the time range of run 1 is known for
entries changed then).

| name | CID | entry | old | new | new / old | action | time |
|---|---|---|---|---|---|---|---|
| (+)-neomenthol | 439263 | value | 0.0707 |  |  | delete | 2026-10-23 12:16-12:27 |
| (1R,2R,4R)-1,7,7-trimethylbicyclo[2.2.1]heptan-2-ol | 6321405 | value | 0.00435 | 0.004 | 0.919 | update | 2026-10-23 12:16-12:27 |
| (E)-2-butenoic acid | 637090 | value | 0.024 | 0.0244 | 1.02 | update | 2026-10-23 12:16-12:27 |
| (E)-2-undecene | 5364452 | null |  | 0.699 |  | update | 2026-10-23 12:40:28 |
| (E)-hex-2-ene | 639661 | value | 21.7 | 23.1 | 1.06 | update | 2026-10-23 12:16-12:27 |
| (E,E)-2,4-hexadienoic acid | 643460 | value | 0.441 | 4.11e-05 | 9.3e-05 | update | 2026-10-23 12:16-12:27 |
| (E,Z)-farnesol | 1549109 | value | 5.25e-06 |  |  | delete | 2026-10-23 12:16-12:27 |
| (tert-butyl)benzene | 7366 | value | 0.449 | 0.293 | 0.653 | update | 2026-10-23 12:16-12:27 |
| (Z)-1,2-dichloroethylene | 643833 | value | 23.9 | 26.7 | 1.12 | update | 2026-10-23 12:16-12:27 |
| (Z)-2-heptene | 643836 | value | 6.45 |  |  | delete | 2026-10-23 12:16-12:27 |
| (Z)-3-hexenol; (Z)-3-hexen-1-ol; cis-3-hexen-1-ol; 3-hexen-1-ol, (Z)- | 5281167 | manual value | 0.125 |  |  | delete | 2026-10-23 12:52:33 |
| (Z)-9-octadecen-1-ol | 5284499 | value | 20.8 | 1.24e-05 | 5.95e-07 | update | 2026-10-23 12:16-12:27 |
| (Z)-farnesol; (Z,E)-farnesol | 1549108 | value | 5.25e-06 |  |  | delete | 2026-10-23 12:16-12:27 |
| (Z,Z)-farnesol | 1549107 | value | 5.25e-06 |  |  | delete | 2026-10-23 12:16-12:27 |
| 1,1,2,2-tetrachloroethane | 6591 | value | 0.718 | 0.765 | 1.07 | update | 2026-10-23 12:16-12:27 |
| 1,1,2-trichloroethane | 6574 | value | 2.73 | 3.07 | 1.12 | update | 2026-10-23 12:16-12:27 |
| 1,1,2-trichloroethene | 6575 | value | 8.65 | 9.2 | 1.06 | update | 2026-10-23 12:16-12:27 |
| 1,1-dichloroethane | 6365 | value | 24.9 | 30.3 | 1.22 | update | 2026-10-23 12:16-12:27 |
| 1,1-diethoxyethane; 1, 1-diethoxyethane; ethane, 1,1-diethoxy- | 7765 | value | 4.41 | 3.68 | 0.834 | update | 2026-10-23 12:16-12:27 |
| 1,1-dimethoxyethane | 10795 | value | 22.8 | 22.8 | 1 | update | 2026-10-23 12:16-12:27 |
| 1,2,3,4-tetrahydronaphthalene | 8404 | value | 0.0711 | 0.0495 | 0.697 | update | 2026-10-23 12:16-12:27 |
| 1,2,3-benzenetriol | 1057 | value | 0.888 | 6.39e-05 | 7.19e-05 | update | 2026-10-23 12:16-12:27 |
| 1,2,4,5-tetrachlorobenzene | 7270 | value | 1.33 | 0.000713 | 0.000535 | update | 2026-10-23 12:16-12:27 |
| 1,2,4-trichlorobenzene | 13 | value | 0.0938 | 0.0542 | 0.578 | update | 2026-10-23 12:16-12:27 |
| 1,2,4-trimethylbenzene | 7247 | value | 0.337 | 0.28 | 0.832 | update | 2026-10-23 12:16-12:27 |
| 1,2-diaminopropane | 6567 | value | 1.2 | 1.2 | 1 | update | 2026-10-23 12:16-12:27 |
| 1,2-dibromoethane | 7839 | value | 1.5 | 1.48 | 0.984 | update | 2026-10-23 12:16-12:27 |
| 1,2-dichlorobenzene | 7239 | value | 0.171 | 0.181 | 1.06 | update | 2026-10-23 12:16-12:27 |
| 1,2-dichloroethane | 11 | value | 10 | 10.5 | 1.05 | update | 2026-10-23 12:16-12:27 |
| 1,2-diethoxyethane | 12375 | value | 2.06 | 0.449 | 0.218 | update | 2026-10-23 12:16-12:27 |
| 1,2-dimethylhydrazine | 1322 | value | 9.21 | 9.16 | 0.994 | update | 2026-10-23 12:16-12:27 |
| 1,2-ethanediamine | 3301 | value | 1.5 | 1.61 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 1,2-ethanediol | 174 | value | 0.00829 | 0.0123 | 1.48 | update | 2026-10-23 12:16-12:27 |
| 1,2-ethanediol diacetate | 8121 | value | 0.0316 | 0.0103 | 0.327 | update | 2026-10-23 12:16-12:27 |
| 1,2-ethanedithiol | 10902 | value | 0.748 | 0.753 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 1,2-propanediol | 1030 | value | 0.0147 | 0.0173 | 1.18 | update | 2026-10-23 12:16-12:27 |
| 1,2-propanediol diacetate | 12198 | value | 0.053 | 0.03 | 0.566 | update | 2026-10-23 12:16-12:27 |
| 1,3,3-trimethyl-2-[(1E,3E,5E,7E,9E,11E,13E,15E,17E)-3,7,12,16-tetramethyl-18-(2,6,6-trimethylcyclohexen-1-yl)octadeca-1,3,5,7,9,11,13,15,17-nonaenyl]cyclohexene | 5280489 | value | 0.024 |  |  | delete | 2026-10-23 12:16-12:27 |
| 1,3,5-trimethylbenzene | 7947 | value | 0.282 | 0.331 | 1.17 | update | 2026-10-23 12:16-12:27 |
| 1,3-benzenediol | 5054 | value | 0.0191 | 3.95e-05 | 0.00207 | update | 2026-10-23 12:16-12:27 |
| 1,3-butanediol | 7896 | value | 0.00444 | 0.00622 | 1.4 | update | 2026-10-23 12:16-12:27 |
| 1,3-dichlorobenzene | 10943 | value | 0.248 | 0.286 | 1.15 | update | 2026-10-23 12:16-12:27 |
| 1,3-diethylbenzene | 8864 | value | 0.155 | 0.16 | 1.03 | update | 2026-10-23 12:16-12:27 |
| 1,3-dimethylbenzene | 7929 | value | 1.11 | 1.11 | 0.997 | update | 2026-10-23 12:16-12:27 |
| 1,3-propanediol | 10442 | value | 0.00561 | 0.00588 | 1.05 | update | 2026-10-23 12:16-12:27 |
| 1,4-benzenedicarbaldehyde | 12173 | value | 0.000536 |  |  | delete | 2026-10-23 12:16-12:27 |
| 1,4-benzenediol | 785 | value | 0.178 | 2.53e-06 | 1.42e-05 | update | 2026-10-23 12:16-12:27 |
| 1,4-cineole | 10106 | value | 0.257 |  |  | delete | 2026-10-23 12:16-12:27 |
| 1,4-diethylbenzene | 7734 | value | 0.139 | 0.137 | 0.986 | update | 2026-10-23 12:16-12:27 |
| 1,4-dimethoxybenzene | 9016 | value | 0.0111 | 0.0107 | 0.958 | update | 2026-10-23 12:16-12:27 |
| 1,4-dimethylbenzene | 7809 | value | 1.16 | 1.18 | 1.02 | update | 2026-10-23 12:16-12:27 |
| 1,4-hexadiene | 5365552 | value | 23 |  |  | delete | 2026-10-23 12:16-12:27 |
| 1,4-naphthalenedione | 8530 | value | 0.00131 | 2.4e-05 | 0.0183 | update | 2026-10-23 12:16-12:27 |
| 1,5-pentanedial | 3485 | value | 1.54 | 1.19 | 0.771 | update | 2026-10-23 12:16-12:27 |
| 1,6-diaminohexane | 16402 | value | 0.133 | 0.0533 | 0.4 | update | 2026-10-23 12:16-12:27 |
| 1,6-dimethylnaphthalene | 11328 | value | 0.00164 | 0.00195 | 1.19 | update | 2026-10-23 12:16-12:27 |
| 1-allyl-4-methoxybenzene | 8815 | value | 0.0143 | 0.00667 | 0.465 | update | 2026-10-23 12:16-12:27 |
| 1-amino-2-propanol | 4 | value | 0.203 | 0.0627 | 0.308 | update | 2026-10-23 12:16-12:27 |
| 1-bromo-2-chloroethane | 7849 | value | 4.72 | 4.41 | 0.935 | update | 2026-10-23 12:16-12:27 |
| 1-bromo-2-phenylethane | 7666 | value | 0.0323 | 0.0326 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 1-bromo-2-propanone | 11715 | value | 4.77 | 1.17 | 0.245 | update | 2026-10-23 12:16-12:27 |
| 1-bromoethane | 6332 | value | 54.8 | 62.3 | 1.14 | update | 2026-10-23 12:16-12:27 |
| 1-bromopropane | 7840 | value | 23.2 | 14.8 | 0.637 | update | 2026-10-23 12:16-12:27 |
| 1-butanamine | 8007 | value | 11 | 12.4 | 1.13 | update | 2026-10-23 12:16-12:27 |
| 1-butanethiol | 8012 | value | 5.02 | 6.07 | 1.21 | update | 2026-10-23 12:16-12:27 |
| 1-butanol; butyl alcohol | 263 | manual value | 1.23 | 0.933 | 0.759 | update | 2026-10-23 12:49:41 |
| 1-chloro-2-methylbenzene | 7238 | value | 0.611 | 0.508 | 0.832 | update | 2026-10-23 12:16-12:27 |
| 1-decanol | 8174 | value | 0.00109 | 0.00113 | 1.04 | update | 2026-10-23 12:16-12:27 |
| 1-decene | 13381 | value | 0.225 | 0.223 | 0.989 | update | 2026-10-23 12:16-12:27 |
| 1-ethoxy-2-propanol | 15287 | value | 1.15 | 1 | 0.867 | update | 2026-10-23 12:16-12:27 |
| 1-heptanethiol | 15422 | value | 0.172 | 0.17 | 0.987 | update | 2026-10-23 12:16-12:27 |
| 1-hexanamine | 8102 | value | 1 | 0.934 | 0.934 | update | 2026-10-23 12:16-12:27 |
| 1-hexanethiol | 8106 | value | 0.16 |  |  | delete | 2026-10-23 12:16-12:27 |
| 1-hydroxy-2-propanone; 2-propanone, 1-hydroxy- | 8299 | value | 0.393 | 0.293 | 0.746 | update | 2026-10-23 12:16-12:27 |
| 1-methoxy-2-propanol | 7900 | value | 1.55 | 1.64 | 1.06 | update | 2026-10-23 12:16-12:27 |
| 1-methyl-2-pyrrolidone | 13387 | value | 0.0434 | 0.0437 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 1-nonanol | 8914 | value | 0.0152 | 0.00303 | 0.199 | update | 2026-10-23 12:16-12:27 |
| 1-octadecanol | 8221 | value | 0.0366 | 3.47e-07 | 9.48e-06 | update | 2026-10-23 12:16-12:27 |
| 1-octanethiol | 8144 | value | 0.195 | 0.0577 | 0.297 | update | 2026-10-23 12:16-12:27 |
| 1-pentanethiol | 8067 | value | 1.85 | 1.85 | 1 | update | 2026-10-23 12:16-12:27 |
| 1-phenyl-2-propanone | 7678 | value | 0.0253 | 0.0293 | 1.16 | update | 2026-10-23 12:16-12:27 |
| 1-phenylethanol | 7409 | value | 0.01 | 0.00667 | 0.667 | update | 2026-10-23 12:16-12:27 |
| 1-propanamine | 7852 | value | 38.9 | 41.3 | 1.06 | update | 2026-10-23 12:16-12:27 |
| 1-propanol | 1031 | manual value | 2.93 | 2.8 | 0.955 | update | 2026-10-23 12:49:35 |
| 1-propene | 8252 | value | 1.11e+03 | 1.16e+03 | 1.04 | update | 2026-10-23 12:16-12:27 |
| 1-tridecanol | 8207 | value | 0.161 | 5.81e-05 | 0.000361 | update | 2026-10-23 12:16-12:27 |
| 1-tridecene | 17095 | value | 0.00825 | 0.00851 | 1.03 | update | 2026-10-23 12:16-12:27 |
| 2(5h)-furanone | 10341 | value | 0.131 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,2'-oxybisethanol | 8117 | value | 0.0278 | 0.00076 | 0.0274 | update | 2026-10-23 12:16-12:27 |
| 2,2,4-trimethylpentane | 10907 | value | 5.84 | 5.26 | 0.901 | update | 2026-10-23 12:16-12:27 |
| 2,2-dimethylpropane | 10041 | value | 159 | 172 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 2,2-dimethylpropanoic acid | 6417 | value | 0.131 | 0.0667 | 0.507 | update | 2026-10-23 12:16-12:27 |
| 2,3,6-trimethylphenol | 17016 | value | 0.00433 | 0.00467 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 2,3-butanediol | 262 | manual value | 0.0322 | 0.0324 | 1.01 | update | 2026-10-23 12:49:43 |
| 2,3-dimethylbutane | 6589 | value | 29 | 31.3 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 2,4,4,6,6-pentamethylhept-2-ene | 170244 | value | 0.52 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,4,6-tribromophenol | 1483 | value | 2.93e-05 | 7e-06 | 0.239 | update | 2026-10-23 12:16-12:27 |
| 2,4,6-trimethyl-1,3,5-trioxane | 31264 | value | 1.95 | 2.19 | 1.12 | update | 2026-10-23 12:16-12:27 |
| 2,4-dibromophenol | 12005 | value | 0.000179 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,4-dimethylphenol | 7771 | value | 0.0108 | 0.0136 | 1.26 | update | 2026-10-23 12:16-12:27 |
| 2,4-dimethylpyridine | 7936 | value | 0.423 | 0.753 | 1.78 | update | 2026-10-23 12:16-12:27 |
| 2,4-hexadien-1-ol | 641256 | value | 0.00267 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,4-pentanedione | 31261 | value | 0.573 | 0.93 | 1.62 | update | 2026-10-23 12:16-12:27 |
| 2,5-dimethylfuran; furan, 2,5-dimethyl- | 12266 | value | 5.26 | 3.45 | 0.657 | update | 2026-10-23 12:16-12:27 |
| 2,5-dimethylphenol | 7267 | value | 0.058 | 0.0208 | 0.358 | update | 2026-10-23 12:16-12:27 |
| 2,5-hexanedione | 8035 | value | 0.302 | 0.391 | 1.29 | update | 2026-10-23 12:16-12:27 |
| 2,6,10,14-tetramethylhexadecane | 12523 | value | 0.00044 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,6,10,14-tetramethylpentadecane | 15979 | value | 0.000587 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,6-dibromophenol | 11847 | value | 0.000179 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2,6-dimethyl-2,5-heptadien-4-one | 10438 | value | 0.0508 | 0.051 | 1 | update | 2026-10-23 12:16-12:27 |
| 2,6-dimethyl-4-heptanol | 7957 | value | 0.227 | 0.04 | 0.176 | update | 2026-10-23 12:16-12:27 |
| 2,6-dimethyl-4-heptanone | 7958 | value | 0.254 | 0.22 | 0.867 | update | 2026-10-23 12:16-12:27 |
| 2,6-dimethylnaphthalene | 11387 | null |  | 0.00112 |  | update | 2026-10-23 12:41:00 |
| 2,6-dimethylphenol | 11335 | value | 0.0363 | 0.0365 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 2-(2-ethoxyethoxy)ethyl acetate | 8165 | value | 0.0126 | 0.0132 | 1.05 | update | 2026-10-23 12:16-12:27 |
| 2-(dimethylamino)ethanol | 7902 | value | 0.625 | 0.424 | 0.678 | update | 2026-10-23 12:16-12:27 |
| 2-[2-(2-ethoxyethoxy)ethoxy]ethanol | 8190 | value | 0.000299 | 0.0003 | 1 | update | 2026-10-23 12:16-12:27 |
| 2-[2-(2-hydroxyethoxy)ethoxy]ethanol | 8172 | value | 0.0334 | 2e-05 | 0.000598 | update | 2026-10-23 12:16-12:27 |
| 2-amino-3,4-dimethylimidazo(4,5-f)quinoline | 62274 | value | 0.132 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2-aminobenzoic acid | 227 | value | 5.81e-05 | 3.71e-05 | 0.638 | update | 2026-10-23 12:16-12:27 |
| 2-aminoethanol | 700 | value | 0.16 | 0.0539 | 0.337 | update | 2026-10-23 12:16-12:27 |
| 2-aminopyridine | 10439 | value | 0.226 | 0.229 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 2-benzylideneheptanal | 31209 | value | 0.000575 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2-bromophenol | 7244 | value | 0.00449 | 0.004 | 0.892 | update | 2026-10-23 12:16-12:27 |
| 2-butanamine | 24874 | value | 26.9 | 23.7 | 0.883 | update | 2026-10-23 12:16-12:27 |
| 2-butanethiol | 10560 | value | 10.8 | 10.8 | 1 | update | 2026-10-23 12:16-12:27 |
| 2-butanol | 6568 | value | 1.96 | 2.44 | 1.25 | update | 2026-10-23 12:16-12:27 |
| 2-butenal; (E)-2-butenal | 447466 | value | 3.27 | 4.36 | 1.33 | update | 2026-10-23 12:16-12:27 |
| 2-butene | 62695 | value | 212 | 227 | 1.07 | update | 2026-10-23 12:16-12:27 |
| 2-butoxyethanol | 8133 | value | 0.111 | 0.112 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 2-butoxyethyl acetate | 8160 | value | 0.0417 | 0.0437 | 1.05 | update | 2026-10-23 12:16-12:27 |
| 2-chlorobiphenyl | 249266 | null |  | 0.000541 |  | update | 2026-10-23 12:42:46 |
| 2-chloronaphthalene | 7056 | value | 0.00156 | 0.00142 | 0.911 | update | 2026-10-23 12:16-12:27 |
| 2-chloropropane | 6361 | value | 68.7 | 68.7 | 1 | update | 2026-10-23 12:16-12:27 |
| 2-ethoxyethanol | 8076 | value | 0.582 | 0.708 | 1.22 | update | 2026-10-23 12:16-12:27 |
| 2-ethyl-1-butanol; 1-butanol, 2-ethyl- | 7358 | value | 0.343 | 0.204 | 0.594 | update | 2026-10-23 12:16-12:27 |
| 2-ethyl-1-hexanol; 1-hexanol, 2-ethyl- | 7720 | value | 0.0194 | 0.0181 | 0.937 | update | 2026-10-23 12:16-12:27 |
| 2-ethylhexanal | 31241 | value | 0.237 | 0.267 | 1.12 | update | 2026-10-23 12:16-12:27 |
| 2-ethylhexanoic acid | 8697 | value | 0.004 | 0.004 | 1 | update | 2026-10-23 12:16-12:27 |
| 2-heptanol | 10976 | value | 0.148 | 0.133 | 0.896 | update | 2026-10-23 12:16-12:27 |
| 2-hexanone | 11583 | value | 1.15 | 1.55 | 1.35 | update | 2026-10-23 12:16-12:27 |
| 2-hydroxybenzoic acid | 338 | value | 0.038 | 1.55e-05 | 0.000407 | update | 2026-10-23 12:16-12:27 |
| 2-isopropoxyethanol; ethanol, 2-(1-methylethoxy)- | 7996 | value | 0.363 | 0.44 | 1.21 | update | 2026-10-23 12:16-12:27 |
| 2-mercaptoethanol | 1567 | value | 0.158 | 0.234 | 1.48 | update | 2026-10-23 12:16-12:27 |
| 2-methoxyethanol | 8019 | value | 1.05 | 1.27 | 1.21 | update | 2026-10-23 12:16-12:27 |
| 2-methoxyethyl acetate | 8054 | value | 0.515 | 0.434 | 0.843 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-1-butanol; 1-butanol, 2-methyl- | 8723 | value | 0.418 | 0.419 | 1 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-1-pentanol; 1-pentanol, 2-methyl- | 7745 | value | 0.22 | 0.256 | 1.17 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-1-propanamine | 6558 | value | 15.9 | 18.4 | 1.16 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-1-propene | 8255 | value | 282 | 308 | 1.09 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2,4-pentanediol | 7870 | value | 0.00673 | 0.00553 | 0.822 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2-butanol | 6405 | value | 1.89 | 2.23 | 1.18 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2-propanamine | 6385 | value | 58.3 | 49.6 | 0.85 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2-propanethiol | 6387 | value | 22.4 | 24.1 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2-propanol | 6386 | value | 5.29 | 5.49 | 1.04 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2-propenal | 6562 | value | 19.1 | 20.7 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-2-propenoic acid | 4093 | value | 0.115 | 0.132 | 1.14 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-3-butyn-2-ol | 8258 | value | 2.07 | 2 | 0.968 | update | 2026-10-23 12:16-12:27 |
| 2-methyl-5-nitroaniline | 7444 | value | 0.00013 | 1.73e-06 | 0.0133 | update | 2026-10-23 12:16-12:27 |
| 2-methylaniline | 7242 | value | 0.0325 | 0.0346 | 1.06 | update | 2026-10-23 12:16-12:27 |
| 2-methylbenzaldehyde | 10722 | value | 0.074 | 0.044 | 0.595 | update | 2026-10-23 12:16-12:27 |
| 2-methylbiphenyl | 12563 | value | 0.0213 | 0.04 | 1.88 | update | 2026-10-23 12:16-12:27 |
| 2-methylbutane | 6556 | value | 85.5 | 79.2 | 0.926 | update | 2026-10-23 12:16-12:27 |
| 2-methyldecane | 23415 | null |  | 0.0653 |  | update | 2026-10-23 12:40:50 |
| 2-methylheptane | 11594 | value | 4.03 | 2.76 | 0.685 | update | 2026-10-23 12:16-12:27 |
| 2-methylhexane; hexane, 2-methyl- | 11582 | value | 7.05 | 8.8 | 1.25 | update | 2026-10-23 12:16-12:27 |
| 2-methylnaphthalene | 7055 | value | 0.007 | 0.00733 | 1.05 | update | 2026-10-23 12:16-12:27 |
| 2-methylpentane | 7892 | value | 27.9 | 28.1 | 1.01 | update | 2026-10-23 12:16-12:27 |
| 2-methylpropanal; propanal, 2-methyl- | 6561 | value | 21 | 23.1 | 1.1 | update | 2026-10-23 12:16-12:27 |
| 2-methylpropanenitrile | 6559 | value | 9.74 | 4.36 | 0.447 | update | 2026-10-23 12:16-12:27 |
| 2-methylpropanol; isobutyl alcohol; 2-methyl-1-propanol; 1-propanol, 2-methyl- | 6560 | value | 1.28 | 1.39 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 2-naphthol | 8663 | value | 0.441 | 0.00102 | 0.00232 | update | 2026-10-23 12:16-12:27 |
| 2-nonylphenol | 67296 | value | 7.87e-06 | 1.25e-05 | 1.59 | update | 2026-10-23 12:16-12:27 |
| 2-octanol | 20083 | value | 0.0574 | 0.0322 | 0.561 | update | 2026-10-23 12:16-12:27 |
| 2-pentanol; methylbutanol | 22386 | manual value | 0.499 | 0.55 | 1.1 | update | 2026-10-23 12:51:01 |
| 2-pentanone; propyl methyl ketone | 7895 | manual value | 3.06 | 4.72 | 1.54 | update | 2026-10-23 12:50:28 |
| 2-phenoxyethanol | 31236 | value | 0.0341 | 0.00247 | 0.0723 | update | 2026-10-23 12:16-12:27 |
| 2-phenyl-2-propanol | 12053 | value | 0.0693 |  |  | delete | 2026-10-23 12:16-12:27 |
| 2-phenylphenol | 7017 | value | 0.0335 | 0.000267 | 0.00796 | update | 2026-10-23 12:16-12:27 |
| 2-propanamine | 6363 | value | 77.7 | 77.3 | 0.994 | update | 2026-10-23 12:16-12:27 |
| 2-propanethiol | 6364 | value | 37 | 37 | 1 | update | 2026-10-23 12:16-12:27 |
| 2-propanol | 3776 | value | 5.11 | 6.05 | 1.19 | update | 2026-10-23 12:16-12:27 |
| 2-propen-1-ol | 7858 | value | 2.83 | 3.28 | 1.16 | update | 2026-10-23 12:16-12:27 |
| 2-propenal | 7847 | manual value | 31.2 | 36.5 | 1.17 | update | 2026-10-23 12:49:22 |
| 2-propenenitrile | 7855 | value | 12.6 | 14.5 | 1.15 | update | 2026-10-23 12:16-12:27 |
| 2-propenoic acid | 6581 | value | 0.454 | 0.529 | 1.17 | update | 2026-10-23 12:16-12:27 |
| 2-tert-butyl-4-methoxyphenol | 8456 | null |  | 0.000331 |  | update | 2026-10-23 12:39:34 |
| 3,4-dimethylphenol | 7249 | value | 0.0474 | 0.00475 | 0.1 | update | 2026-10-23 12:16-12:27 |
| 3,5,5-trimethyl-1-hexanol | 18938 | value | 0.035 | 0.03 | 0.857 | update | 2026-10-23 12:16-12:27 |
| 3,5,5-trimethyl-2-cyclohexen-1-one | 6544 | value | 0.0437 | 0.0584 | 1.34 | update | 2026-10-23 12:16-12:27 |
| 3,5-dimethylphenol | 7948 | value | 0.0373 | 0.0052 | 0.14 | update | 2026-10-23 12:16-12:27 |
| 3-(methylthio)propanal | 18635 | value | 0.0843 | 0.0529 | 0.628 | update | 2026-10-23 12:16-12:27 |
| 3-buten-2-one | 6570 | value | 13.1 | 15.6 | 1.2 | update | 2026-10-23 12:16-12:27 |
| 3-cyclohexene-1-methanol, .alpha.,.alpha.,4-trimethyl-, (s)- | 443162 | null |  | 0.00667 |  | update | 2026-10-23 12:37:26 |
| 3-heptanone | 7802 | value | 0.413 | 0.187 | 0.452 | update | 2026-10-23 12:16-12:27 |
| 3-hydroxy-2,2,4-trimethylpentyl 2-methylpropanoate | 6490 | value | 0.0013 | 0.0013 | 1 | update | 2026-10-23 12:16-12:27 |
| 3-hydroxybutanal | 7897 | value | 1.92 | 2.8 | 1.45 | update | 2026-10-23 12:16-12:27 |
| 3-methyl-2-butanone; 3-methylbutan-2-one | 11251 | value | 6.52 | 5.5 | 0.843 | update | 2026-10-23 12:16-12:27 |
| 3-methyl-2-buten-1-ol | 11173 | manual value | 0.14 | 0.14 | 1 | update | 2026-10-23 12:50:31 |
| 3-methyl-2-butenal | 61020 | value | 0.748 | 0.75 | 1 | update | 2026-10-23 12:16-12:27 |
| 3-methylbenzaldehyde | 12105 | value | 0.0505 |  |  | delete | 2026-10-23 12:16-12:27 |
| 3-methylbenzoic acid | 7418 | value | 0.0667 | 1.87e-05 | 0.00028 | update | 2026-10-23 12:16-12:27 |
| 3-methylbutanal; butanal, 3-methyl- | 11552 | value | 6.38 | 6.1 | 0.956 | update | 2026-10-23 12:16-12:27 |
| 3-methylbutanol; isoamyl alcohol; iso-amyl alcohol; isopentyl alcohol; 3-methyl-1-butanol; 1-butanol, 3-methyl- | 31260 | manual value | 1.7 | 0.316 | 0.186 | update | 2026-10-23 12:51:03 |
| 3-methylpentane | 7282 | value | 17.8 | 25.3 | 1.42 | update | 2026-10-23 12:16-12:27 |
| 3-methylphenol | 342 | value | 0.0133 | 0.0167 | 1.26 | update | 2026-10-23 12:16-12:27 |
| 3-methylpyridine | 7970 | value | 0.738 | 0.807 | 1.09 | update | 2026-10-23 12:16-12:27 |
| 3-pentanol | 11428 | value | 0.985 | 0.8 | 0.812 | update | 2026-10-23 12:16-12:27 |
| 3-pentanone; pentan-3-one | 7288 | value | 4.14 | 4.79 | 1.16 | update | 2026-10-23 12:16-12:27 |
| 3-penten-2-ol; (E)-3-penten-2-ol | 5366239 | value | 0.801 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-(2-furyl)-3-buten-2-one; (E)-4-(2-furyl)-3-buten-2-one | 735940 | value | 0.044 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-(3-methylbut-2-enoxy)furo[3,2-g]chromen-7-one | 68081 | value | 1.87e-08 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-(4-hydroxyphenyl)-2-butanone | 21648 | value | 7.6e-05 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-aminoacetophenone | 7468 | value | 0.000688 | 0.00132 | 1.92 | update | 2026-10-23 12:16-12:27 |
| 4-biphenylamine | 7102 | value | 0.0762 | 5.98e-05 | 0.000785 | update | 2026-10-23 12:16-12:27 |
| 4-bromophenol | 7808 | value | 0.00145 | 0.00156 | 1.08 | update | 2026-10-23 12:16-12:27 |
| 4-chlorophenol | 4684 | value | 0.0121 | 0.0116 | 0.955 | update | 2026-10-23 12:16-12:27 |
| 4-ethylphenol | 31242 | value | 0.00448 | 0.00496 | 1.11 | update | 2026-10-23 12:16-12:27 |
| 4-hydroxy-4-methyl-2-pentanone | 31256 | value | 0.217 | 0.127 | 0.587 | update | 2026-10-23 12:16-12:27 |
| 4-hydroxybutanoic acid | 10413 | value | 0.0173 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-hydroxynonanoic acid lactone; dihydro-5-pentyl-2(3h)-furanone | 7710 | value | 0.00157 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-methoxyphenol | 9015 | value | 0.001 | 0.0009 | 0.897 | update | 2026-10-23 12:16-12:27 |
| 4-methyl-2-pentanone; 2-pentanone, 4-methyl- | 7909 | value | 2.17 | 2.65 | 1.22 | update | 2026-10-23 12:16-12:27 |
| 4-methyl-3-penten-2-one | 8858 | value | 1.16 | 1.09 | 0.94 | update | 2026-10-23 12:16-12:27 |
| 4-methylaniline | 7813 | value | 0.106 | 0.0381 | 0.361 | update | 2026-10-23 12:16-12:27 |
| 4-methylbenzenethiol | 7811 | value | 0.107 | 0.108 | 1 | update | 2026-10-23 12:16-12:27 |
| 4-methylbenzyl alcohol | 11505 | value | 0.0042 | 0.00707 | 1.68 | update | 2026-10-23 12:16-12:27 |
| 4-methyloxetan-2-one | 18303 | value | 0.143 | 0.22 | 1.53 | update | 2026-10-23 12:16-12:27 |
| 4-methylpyridine | 7963 | value | 0.766 | 0.765 | 0.998 | update | 2026-10-23 12:16-12:27 |
| 4-nonylphenol | 1752 | value | 0.000109 | 0.000109 | 1 | update | 2026-10-23 12:16-12:27 |
| 4-terpineol; terpinen-4-ol | 11230 | value | 0.00533 |  |  | delete | 2026-10-23 12:16-12:27 |
| 4-tert-butylphenol | 7393 | value | 0.133 | 0.00123 | 0.00919 | update | 2026-10-23 12:16-12:27 |
| 4-vinylpyridine | 7502 | value | 0.241 | 0.247 | 1.03 | update | 2026-10-23 12:16-12:27 |
| 5,8-dimethoxypsoralen | 68079 | value | 9.07e-08 |  |  | delete | 2026-10-23 12:16-12:27 |
| 5-(1-propenyl)-1,3-benzodioxole | 637796 | value | 0.0314 | 0.00267 | 0.0848 | update | 2026-10-23 12:16-12:27 |
| 5-allyl-1,3-benzodioxole | 5144 | value | 0.071 | 0.00941 | 0.133 | update | 2026-10-23 12:16-12:27 |
| 5-ethyl-2-methylpyridine | 7728 | value | 0.315 | 0.191 | 0.604 | update | 2026-10-23 12:16-12:27 |
| 5-methoxypsoralen | 2355 | value | 1e-07 |  |  | delete | 2026-10-23 12:16-12:27 |
| 5-methyl-2-hexanone | 8034 | value | 0.69 | 0.769 | 1.12 | update | 2026-10-23 12:16-12:27 |
| 5-methyl-3-heptanone | 7822 | value | 0.267 | 0.267 | 1 | update | 2026-10-23 12:16-12:27 |
| 5-nonanone | 10405 | value | 0.0735 | 0.0736 | 1 | update | 2026-10-23 12:16-12:27 |
| 6-hydroxyhexanoic acid lactone | 10401 | value | 0.012 | 0.0094 | 0.78 | update | 2026-10-23 12:16-12:27 |
| 6-methyl-1-heptanol | 15450 | value | 0.0679 | 0.0507 | 0.747 | update | 2026-10-23 12:16-12:27 |
| 6-methylheptyl 2-propenoate | 165746 | value | 0.0766 | 0.133 | 1.74 | update | 2026-10-23 12:16-12:27 |
| 7-methoxy-8-(3-methyl-2-butenyl)coumarin | 10228 | value | 9.2e-07 |  |  | delete | 2026-10-23 12:16-12:27 |
| 8-hydroxyquinoline | 1923 | value | 0.000217 | 0.000221 | 1.02 | update | 2026-10-23 12:16-12:27 |
| 8-methoxypsoralen | 4114 | value | 2.87e-07 |  |  | delete | 2026-10-23 12:16-12:27 |
| 8-methylnonanal | 14862 | value | 0.123 | 0.0387 | 0.315 | update | 2026-10-23 12:16-12:27 |
| 9,12,15-octadecatrienoic acid; (Z,Z,Z)-9,12,15-octadecatrienoic acid | 5280934 | value | 0.00333 | 7.2e-08 | 2.16e-05 | update | 2026-10-23 12:16-12:27 |
| [(E)-prop-1-enyl]benzene | 252325 | value | 0.14 | 0.15 | 1.07 | update | 2026-10-23 12:16-12:27 |
| [2,2,4-trimethyl-3-(2-methylpropanoyloxy)pentyl] 2-methylpropanoate | 23284 | value | 0.00113 |  |  | delete | 2026-10-23 12:16-12:27 |
| acenaphthene | 6734 | value | 0.000293 | 0.000269 | 0.917 | update | 2026-10-23 12:16-12:27 |
| acenaphthylene | 9161 | value | 0.000765 | 0.00064 | 0.836 | update | 2026-10-23 12:16-12:27 |
| acetaldehyde; acetylaldehyde | 177 | manual value | 107 | 120 | 1.13 | update | 2026-10-23 12:49:05 |
| acetamide | 178 | value | 0.101 | 0.00533 | 0.0527 | update | 2026-10-23 12:16-12:27 |
| acetic acid | 176 | value | 1.72 | 2.09 | 1.21 | update | 2026-10-23 12:16-12:27 |
| acetic anhydride; acetic acid, anhydride | 7918 | value | 0.71 | 0.68 | 0.958 | update | 2026-10-23 12:16-12:27 |
| acetone; 2-propanone | 180 | manual value | 26.3 | 30.8 | 1.17 | update | 2026-10-23 12:49:27 |
| acetophenone; ethanone, 1-phenyl- | 7410 | manual value | 0.0793 | 0.0529 | 0.667 | update | 2026-10-23 13:01:36 |
| adenosine | 60961 | value | 0.08 |  |  | delete | 2026-10-23 12:16-12:27 |
| allyl acetate | 11584 | value | 3.63 | 3.68 | 1.01 | update | 2026-10-23 12:16-12:27 |
| allyl isothiocyanate | 5971 | value | 0.696 | 0.538 | 0.773 | update | 2026-10-23 12:16-12:27 |
| allyl propyl disulfide | 16591 | value | 0.0513 | 0.052 | 1.01 | update | 2026-10-23 12:16-12:27 |
| allylamine | 7853 | value | 30.3 | 29.5 | 0.973 | update | 2026-10-23 12:16-12:27 |
| alpha-pinene; .alpha.-pinene | 6654 | value | 0.867 | 0.633 | 0.731 | update | 2026-10-23 12:16-12:27 |
| ammonia | 222 | value | 934 | 1.01e+03 | 1.08 | update | 2026-10-23 12:16-12:27 |
| aniline | 6115 | value | 0.0772 | 0.0891 | 1.15 | update | 2026-10-23 12:16-12:27 |
| anthracene | 8418 | value | 0.0667 | 4.04e-05 | 0.000606 | update | 2026-10-23 12:16-12:27 |
| apigenin | 5280443 | value | 0.0135 |  |  | delete | 2026-10-23 12:16-12:27 |
| ascorbic acid | 54670067 | value | 0.124 |  |  | delete | 2026-10-23 12:16-12:27 |
| azepan-2-one | 7768 | value | 0.114 | 0.000257 | 0.00224 | update | 2026-10-23 12:16-12:27 |
| benz[a]anthracene | 5954 | value | 0.0973 | 2.8e-08 | 2.88e-07 | update | 2026-10-23 12:16-12:27 |
| benzaldehyde | 240 | manual value | 0.113 | 0.106 | 0.942 | update | 2026-10-23 12:56:10 |
| benzamide | 2331 | value | 0.000125 | 7.36e-05 | 0.588 | update | 2026-10-23 12:16-12:27 |
| benzene | 241 | value | 10.9 | 12.6 | 1.16 | update | 2026-10-23 12:16-12:27 |
| benzenediol; 1,2-benzenediol | 289 | value | 0.779 | 0.003 | 0.00385 | update | 2026-10-23 12:16-12:27 |
| benzeneethanol; 2-phenylethanol; phenethyl alcohol; phenylethyl alcohol | 6054 | manual value | 0.0101 | 0.0116 | 1.15 | update | 2026-10-23 12:56:52 |
| benzenethiol | 7969 | value | 0.168 | 0.262 | 1.56 | update | 2026-10-23 12:16-12:27 |
| benzo[a]fluorene | 9195 | value | 3.28e-07 | 3.24e-08 | 0.0987 | update | 2026-10-23 12:16-12:27 |
| benzo[ghi]perylene | 9117 | value | 0.0133 | 0.00667 | 0.5 | update | 2026-10-23 12:16-12:27 |
| benzofuran | 9223 | value | 0.521 | 0.0593 | 0.114 | update | 2026-10-23 12:16-12:27 |
| benzoic acid | 243 | value | 0.0334 | 9.67e-05 | 0.00289 | update | 2026-10-23 12:16-12:27 |
| benzonitrile | 7505 | value | 0.11 | 0.102 | 0.932 | update | 2026-10-23 12:16-12:27 |
| benzoquinone; 1,4-benzoquinone | 4650 | value | 0.013 | 0.0133 | 1.03 | update | 2026-10-23 12:16-12:27 |
| benzyl acetate | 8785 | value | 0.101 | 0.0791 | 0.786 | update | 2026-10-23 12:16-12:27 |
| benzyl alcohol | 244 | manual value | 0.0129 | 0.0125 | 0.973 | update | 2026-10-23 12:56:22 |
| benzyl benzoate | 2345 | value | 2.99e-05 | 2.96e-05 | 0.991 | update | 2026-10-23 12:16-12:27 |
| benzyl chloride | 7503 | value | 0.141 | 0.164 | 1.16 | update | 2026-10-23 12:16-12:27 |
| benzylamine | 7504 | value | 0.0878 | 0.0874 | 0.996 | update | 2026-10-23 12:16-12:27 |
| beta-myrcene; .beta.-myrcene | 31253 | value | 0.387 | 0.279 | 0.72 | update | 2026-10-23 12:16-12:27 |
| betaine | 247 | value | 3.62e-08 | 7.07e-08 | 1.95 | update | 2026-10-23 12:16-12:27 |
| biphenyl | 7095 | value | 0.000928 | 0.00119 | 1.28 | update | 2026-10-23 12:16-12:27 |
| bromochloromethane | 6333 | value | 16.5 | 18.9 | 1.15 | update | 2026-10-23 12:16-12:27 |
| bromomethane | 6323 | value | 763 | 216 | 0.283 | update | 2026-10-23 12:16-12:27 |
| butanal; butyraldehyde | 261 | manual value | 13.5 | 14.8 | 1.1 | update | 2026-10-23 12:50:07 |
| butane | 7843 | value | 223 | 243 | 1.09 | update | 2026-10-23 12:16-12:27 |
| butanenitrile | 8008 | value | 2.19 | 2.6 | 1.19 | update | 2026-10-23 12:16-12:27 |
| butanone; 2-butanone | 6569 | manual value | 11.5 | 12.1 | 1.05 | update | 2026-10-23 12:50:08 |
| butene; 1-butene | 7844 | value | 423 | 300 | 0.71 | update | 2026-10-23 12:16-12:27 |
| butyl 2-methyl-2-propenoate | 7354 | value | 0.276 | 0.283 | 1.02 | update | 2026-10-23 12:16-12:27 |
| butyl acetate; butyl ethanoate; acetic acid, butyl ester | 31272 | manual value | 1.38 | 1.53 | 1.11 | update | 2026-10-23 12:52:40 |
| butyl butyrate; butyl butanoate; butanoic acid, butyl ester | 7983 | value | 0.999 | 0.81 | 0.81 | update | 2026-10-23 12:16-12:27 |
| butyl formate; formic acid, butyl ester | 11614 | value | 3.38 | 2.9 | 0.859 | update | 2026-10-23 12:16-12:27 |
| butyl isocyanate | 8110 | value | 2.26 | 2.1 | 0.927 | update | 2026-10-23 12:16-12:27 |
| butyl propanoate | 11529 | value | 0.485 | 0.38 | 0.784 | update | 2026-10-23 12:16-12:27 |
| butyric acid; butanoic acid | 264 | value | 0.139 | 0.22 | 1.59 | update | 2026-10-23 12:16-12:27 |
| butyrolactone; gamma-butyrolactone; 2(3h)-furanone, dihydro-; 4-hydroxybutanoic acid lactone | 7302 | value | 0.09 | 0.06 | 0.667 | update | 2026-10-23 12:16-12:27 |
| butyrone; 4-heptanone | 31246 | value | 0.681 | 0.696 | 1.02 | update | 2026-10-23 12:16-12:27 |
| caffeine | 2519 | value | 6.07e-08 | 1.2e-07 | 1.98 | update | 2026-10-23 12:16-12:27 |
| camphene | 6616 | value | 0.356 | 0.333 | 0.937 | update | 2026-10-23 12:16-12:27 |
| carbon disulfide | 6348 | value | 46 | 48 | 1.04 | update | 2026-10-23 12:16-12:27 |
| carbondioxide; carbon dioxide | 280 | value | 5.83e+03 | 6.44e+03 | 1.11 | update | 2026-10-23 12:16-12:27 |
| carbonyl sulfide | 10039 | value | 1.25e+03 | 1.25e+03 | 1 | update | 2026-10-23 12:16-12:27 |
| carvone | 7439 | value | 0.0183 | 0.0213 | 1.16 | update | 2026-10-23 12:16-12:27 |
| cedrol | 65575 | value | 1.15e-05 |  |  | delete | 2026-10-23 12:16-12:27 |
| chlorobenzaldehyde | 6996 | value | 0.0987 | 0.04 | 0.405 | update | 2026-10-23 12:16-12:27 |
| chlorobenzene | 7964 | value | 1.39 | 1.6 | 1.15 | update | 2026-10-23 12:16-12:27 |
| chloroethane | 6337 | value | 134 | 134 | 1 | update | 2026-10-23 12:16-12:27 |
| chloromethane | 6327 | value | 540 | 573 | 1.06 | update | 2026-10-23 12:16-12:27 |
| cholesterol | 5997 | value | 0.0228 |  |  | delete | 2026-10-23 12:16-12:27 |
| chrysene | 9171 | value | 4.24e-08 | 8.31e-10 | 0.0196 | update | 2026-10-23 12:16-12:27 |
| chrysophanic acid | 10208 | value | 0.0492 |  |  | delete | 2026-10-23 12:16-12:27 |
| citral; geranial | 638011 | value | 0.0525 | 0.012 | 0.229 | update | 2026-10-23 12:16-12:27 |
| citric acid | 311 | value | 2.27e-09 |  |  | delete | 2026-10-23 12:16-12:27 |
| citronellal | 7794 | value | 0.246 | 0.352 | 1.43 | update | 2026-10-23 12:16-12:27 |
| citronellol; beta-citronellol | 8842 | manual value | 0.004 | 0.00267 | 0.667 | update | 2026-10-23 12:38:20 |
| coumarin | 323 | value | 0.0662 | 0.000131 | 0.00197 | update | 2026-10-23 12:16-12:27 |
| cyclododecane | 9268 | value | 0.00397 | 0.00393 | 0.992 | update | 2026-10-23 12:16-12:27 |
| cyclohexane | 8078 | value | 16.9 | 13 | 0.767 | update | 2026-10-23 12:16-12:27 |
| cyclohexanol | 7966 | value | 0.303 | 0.0876 | 0.289 | update | 2026-10-23 12:16-12:27 |
| cyclohexanone | 7967 | value | 0.643 | 0.583 | 0.907 | update | 2026-10-23 12:16-12:27 |
| cyclohexene | 8079 | value | 9.9 | 11.9 | 1.2 | update | 2026-10-23 12:16-12:27 |
| cyclohexyl acetate | 12146 | value | 1.46 | 1.46 | 0.999 | update | 2026-10-23 12:16-12:27 |
| cyclohexylamine | 7965 | value | 1.42 | 1.35 | 0.951 | update | 2026-10-23 12:16-12:27 |
| cyclohexylmethane; methylcyclohexane; cyclohexane, methyl- | 7962 | value | 6.57 | 5.93 | 0.903 | update | 2026-10-23 12:16-12:27 |
| cyclooctatetraene; 1,3,5,7-cyclooctatetraene | 637866 | value | 1.04 |  |  | delete | 2026-10-23 12:16-12:27 |
| cyclopentadiene | 7612 | value | 55.2 | 45.1 | 0.817 | update | 2026-10-23 12:16-12:27 |
| cyclopentane | 9253 | value | 49 | 42.4 | 0.865 | update | 2026-10-23 12:16-12:27 |
| cyclopentanol | 7298 | value | 126 | 252 | 2 | update | 2026-10-23 12:16-12:27 |
| d,l-methionine | 876 | value | 6.93e-08 |  |  | delete | 2026-10-23 12:16-12:27 |
| d-camphor | 159055 | value | 0.0943 |  |  | delete | 2026-10-23 12:16-12:27 |
| d-limonene; (r)-(+)-limonene; (4r)-(+)-limonene | 440917 | value | 0.213 | 0.264 | 1.24 | update | 2026-10-23 12:16-12:27 |
| decahydronaphthalene; cis-decahydronaphthalene; trans-decahydronaphthalene | 7044 | value | 0.353 | 0.307 | 0.868 | update | 2026-10-23 12:16-12:27 |
| decanal | 8175 | value | 0.0135 | 0.0137 | 1.01 | update | 2026-10-23 12:16-12:27 |
| decane | 15600 | value | 0.228 | 0.184 | 0.807 | update | 2026-10-23 12:16-12:27 |
| decanedioic acid | 5192 | value | 0.0667 | 3.85e-07 | 5.78e-06 | update | 2026-10-23 12:16-12:27 |
| decanoic acid | 2969 | value | 0.0445 | 4.88e-05 | 0.0011 | update | 2026-10-23 12:16-12:27 |
| di-tert-butyl peroxide | 8033 | value | 2.97 | 3.35 | 1.13 | update | 2026-10-23 12:16-12:27 |
| diallylamine | 31279 | value | 2.6 | 2.69 | 1.03 | update | 2026-10-23 12:16-12:27 |
| dibenzo[a,h]anthracene | 5889 | value | 0.0637 | 0.127 | 2 | update | 2026-10-23 12:16-12:27 |
| dibenzofuran | 568 | value | 0.000331 | 0.000416 | 1.26 | update | 2026-10-23 12:16-12:27 |
| dibromo(chloro)methane | 31296 | value | 2.93 | 0.739 | 0.252 | update | 2026-10-23 12:16-12:27 |
| dibromomethane | 3024 | value | 6.72 | 5.92 | 0.881 | update | 2026-10-23 12:16-12:27 |
| dibutyl phthalate | 3026 | value | 0.268 | 2.27e-06 | 8.47e-06 | update | 2026-10-23 12:16-12:27 |
| dibutylamine | 8148 | value | 0.32 | 0.345 | 1.08 | update | 2026-10-23 12:16-12:27 |
| dichlorobenzene; 1,4-dichlorobenzene | 4685 | value | 0.177 | 0.232 | 1.31 | update | 2026-10-23 12:16-12:27 |
| dichloromethane | 6344 | value | 52.6 | 58 | 1.1 | update | 2026-10-23 12:16-12:27 |
| diethyl carbonate | 7766 | value | 1.33 | 1.44 | 1.08 | update | 2026-10-23 12:16-12:27 |
| diethyl ether | 3283 | value | 63.1 | 71.7 | 1.14 | update | 2026-10-23 12:16-12:27 |
| diethyl fumarate | 638144 | value | 0.0599 | 0.0227 | 0.378 | update | 2026-10-23 12:16-12:27 |
| diethyl hexanedioate | 8844 | value | 0.0072 | 0.00773 | 1.07 | update | 2026-10-23 12:16-12:27 |
| diethyl malonate | 7761 | value | 0.03 | 0.0307 | 1.02 | update | 2026-10-23 12:16-12:27 |
| diethyl succinate | 31249 | manual value | 0.0692 | 0.00533 | 0.0771 | update | 2026-10-23 12:59:21 |
| diethylene glycol diethyl ether | 8179 | value | 0.0724 | 0.0693 | 0.957 | update | 2026-10-23 12:16-12:27 |
| diethylene glycol monobutyl ether | 8177 | value | 0.00321 | 0.00462 | 1.44 | update | 2026-10-23 12:16-12:27 |
| diethylene glycol monoethyl ether | 8146 | value | 0.0173 | 0.0179 | 1.04 | update | 2026-10-23 12:16-12:27 |
| diethylene glycol monomethyl ether | 8134 | value | 0.0313 | 0.0333 | 1.06 | update | 2026-10-23 12:16-12:27 |
| diisobutylamine | 8085 | value | 1.06 | 0.97 | 0.914 | update | 2026-10-23 12:16-12:27 |
| dimethoxymethane | 8020 | value | 45.8 | 43.8 | 0.956 | update | 2026-10-23 12:16-12:27 |
| dimethyl carbonate | 12021 | value | 7.39 | 7.39 | 1 | update | 2026-10-23 12:16-12:27 |
| dimethyl disulfide; disulfide, dimethyl | 12232 | value | 3.82 | 3.81 | 0.999 | update | 2026-10-23 12:16-12:27 |
| dimethyl ether | 8254 | value | 396 | 593 | 1.5 | update | 2026-10-23 12:16-12:27 |
| dimethyl sulfide; methane, thiobis- | 1068 | value | 62.4 | 66.9 | 1.07 | update | 2026-10-23 12:16-12:27 |
| dimethyl sulfoxide | 679 | value | 0.0692 | 0.08 | 1.16 | update | 2026-10-23 12:16-12:27 |
| dimethylbenzene; 1,2-dimethylbenzene; benzene, 1,2-dimethyl- | 7237 | value | 0.989 | 0.7 | 0.708 | update | 2026-10-23 12:16-12:27 |
| dimethylnaphthalene; 1,2-dimethylnaphthalene | 11317 | value | 6.4 |  |  | delete | 2026-10-23 12:16-12:27 |
| dimethylphenol; 2,3-dimethylphenol | 10687 | value | 0.372 | 0.0119 | 0.0319 | update | 2026-10-23 12:16-12:27 |
| dioctyl phthalate | 8346 | value | 1.24e-07 | 1.8e-07 | 1.45 | update | 2026-10-23 12:16-12:27 |
| dioxane; 1,4-dioxane | 31275 | value | 4.45 | 5.08 | 1.14 | update | 2026-10-23 12:16-12:27 |
| diphenylamine | 11487 | value | 0.08 | 8.93e-05 | 0.00112 | update | 2026-10-23 12:16-12:27 |
| dl-camphor | 2537 | value | 0.0439 | 0.0867 | 1.97 | update | 2026-10-23 12:16-12:27 |
| dl-menthone | 6986 | value | 0.036 |  |  | delete | 2026-10-23 12:16-12:27 |
| dodecane | 8182 | value | 0.0251 | 0.018 | 0.717 | update | 2026-10-23 12:16-12:27 |
| dodecanoic acid | 3893 | value | 0.133 |  |  | delete | 2026-10-23 12:16-12:27 |
| dodecene; 1-dodecene | 8183 | value | 0.0234 | 0.00206 | 0.0882 | update | 2026-10-23 12:16-12:27 |
| dotriacontane | 11008 | value | 2.27e-08 |  |  | delete | 2026-10-23 12:16-12:27 |
| ellagic acid | 5281855 | value | 0.0375 |  |  | delete | 2026-10-23 12:16-12:27 |
| endo-borneol | 1201518 | value | 0.00535 |  |  | delete | 2026-10-23 12:16-12:27 |
| ethanamine | 6341 | value | 114 | 140 | 1.22 | update | 2026-10-23 12:16-12:27 |
| ethane | 6324 | value | 4.02e+03 | 3.85e+03 | 0.957 | update | 2026-10-23 12:16-12:27 |
| ethanedial | 7860 | value | 2.4 | 34 | 14.2 | update | 2026-10-23 12:16-12:27 |
| ethanethiol | 6343 | value | 62.8 | 70.5 | 1.12 | update | 2026-10-23 12:16-12:27 |
| ethanol | 702 | manual value | 6.67 | 7.91 | 1.19 | update | 2026-10-23 12:49:08 |
| ethene | 6325 | value | 6.36e+03 | 6.89e+03 | 1.08 | update | 2026-10-23 12:16-12:27 |
| ethene, ethoxy-; ethyl vinyl ether | 8023 | value | 64.3 | 68.7 | 1.07 | update | 2026-10-23 12:16-12:27 |
| ethoxybenzene | 7674 | value | 0.206 | 0.208 | 1.01 | update | 2026-10-23 12:16-12:27 |
| ethoxyethyl acetate; 2-ethoxyethyl acetate | 8095 | value | 0.249 | 0.267 | 1.07 | update | 2026-10-23 12:16-12:27 |
| ethyl 2-propenoate | 8821 | value | 4.33 | 5.15 | 1.19 | update | 2026-10-23 12:16-12:27 |
| ethyl 3-oxobutanoate | 8868 | value | 0.102 | 0.104 | 1.02 | update | 2026-10-23 12:16-12:27 |
| ethyl acetate; acetic acid ethyl ester; acetic acid, ethyl ester | 8857 | manual value | 11.3 | 12.4 | 1.1 | update | 2026-10-23 12:50:10 |
| ethyl anthranilate | 6877 | value | 0.000308 | 0.000309 | 1 | update | 2026-10-23 12:16-12:27 |
| ethyl butyrate; ethyl butanoate; butyric acid, ethyl ester; butanoic acid, ethyl ester | 7762 | manual value | 1.79 | 1.87 | 1.04 | update | 2026-10-23 12:52:41 |
| ethyl carbamate | 5641 | value | 0.189 | 0.0415 | 0.22 | update | 2026-10-23 12:16-12:27 |
| ethyl formate; ethyl methanoate; formic acid, ethyl ester | 8025 | manual value | 28.5 | 32.7 | 1.15 | update | 2026-10-23 12:49:28 |
| ethyl lactate; ethyl 2-hydroxypropanoate | 7344 | manual value | 0.409 | 0.5 | 1.22 | update | 2026-10-23 12:50:44 |
| ethyl methacrylate; ethyl 2-methyl-2-propenoate; methacrylic acid, ethyl ester | 7343 | value | 2.76 | 2.5 | 0.905 | update | 2026-10-23 12:16-12:27 |
| ethyl phenylacetate | 7590 | value | 0.00816 | 0.008 | 0.98 | update | 2026-10-23 12:16-12:27 |
| ethyl phthalate; diethyl phthalate | 6781 | value | 0.0269 | 0.000271 | 0.0101 | update | 2026-10-23 12:16-12:27 |
| ethyl propanoate; ethyl propionate; propanoic acid, ethyl ester | 7749 | manual value | 4.96 | 4.77 | 0.962 | update | 2026-10-23 12:50:35 |
| ethyl thiocyanate | 10968 | value | 0.644 | 0.533 | 0.828 | update | 2026-10-23 12:16-12:27 |
| ethylbenzene; benzene, ethyl- | 7500 | value | 1.11 | 1.31 | 1.18 | update | 2026-10-23 12:16-12:27 |
| ethylcyclohexane | 15504 | value | 2.92 | 1.71 | 0.584 | update | 2026-10-23 12:16-12:27 |
| ethylphenol; 2-ethylphenol | 6997 | value | 0.0579 | 0.0204 | 0.352 | update | 2026-10-23 12:16-12:27 |
| eugenol | 3314 | manual value | 0.00276 | 0.00295 | 1.07 | update | 2026-10-23 12:34:41 |
| fluoranthene | 9154 | value | 0.000667 | 1.23e-06 | 0.00184 | update | 2026-10-23 12:16-12:27 |
| fluorene | 6853 | value | 0.444 | 8e-05 | 0.00018 | update | 2026-10-23 12:16-12:27 |
| formaldehyde | 712 | value | 259 | 519 | 2 | update | 2026-10-23 12:16-12:27 |
| formamide | 713 | value | 0.0875 | 0.00813 | 0.093 | update | 2026-10-23 12:16-12:27 |
| formic acid; methanoic acid | 284 | value | 4.99 | 5.68 | 1.14 | update | 2026-10-23 12:16-12:27 |
| furan | 8029 | value | 72.9 | 80 | 1.1 | update | 2026-10-23 12:16-12:27 |
| furfural; 2-furfural | 7362 | manual value | 0.234 | 0.15 | 0.64 | update | 2026-10-23 12:51:19 |
| furfuryl acetate | 12170 | value | 0.133 |  |  | delete | 2026-10-23 12:16-12:27 |
| furfuryl alcohol | 7361 | value | 0.0715 | 0.0867 | 1.21 | update | 2026-10-23 12:16-12:27 |
| genistein | 5280961 | value | 0.0693 |  |  | delete | 2026-10-23 12:16-12:27 |
| geranyl acetate | 1549026 | value | 0.00369 | 0.004 | 1.08 | update | 2026-10-23 12:16-12:27 |
| glycerol | 753 | manual value | 0.000158 | 1.62e-05 | 0.103 | update | 2026-10-23 12:49:36 |
| guaiacol; methoxyphenol; 2-methoxyphenol | 460 | value | 2.23 | 0.0133 | 0.00599 | update | 2026-10-23 12:16-12:27 |
| heptanal | 8130 | value | 0.446 | 0.469 | 1.05 | update | 2026-10-23 12:16-12:27 |
| heptane | 8900 | value | 5.42 | 6.13 | 1.13 | update | 2026-10-23 12:16-12:27 |
| heptanoic acid | 8094 | value | 0.0343 | 0.00109 | 0.0316 | update | 2026-10-23 12:16-12:27 |
| heptanol; 1-heptanol | 8129 | value | 0.0239 | 0.0288 | 1.2 | update | 2026-10-23 12:16-12:27 |
| heptanone; 2-heptanone | 8051 | manual value | 0.405 | 0.357 | 0.88 | update | 2026-10-23 12:55:13 |
| hexachloro-2-propanone | 8303 | value | 0.016 | 0.0501 | 3.13 | update | 2026-10-23 12:16-12:27 |
| hexadecane | 11006 | value | 0.000195 | 0.000199 | 1.02 | update | 2026-10-23 12:16-12:27 |
| hexadecanethiol | 18015 | value | 0.0125 | 0.01 | 0.8 | update | 2026-10-23 12:16-12:27 |
| hexadecanoic acid | 985 | value | 0.0443 | 5.07e-08 | 1.14e-06 | update | 2026-10-23 12:16-12:27 |
| hexadecanol; 1-hexadecanol | 2682 | value | 6.04e-07 | 8e-07 | 1.32 | update | 2026-10-23 12:16-12:27 |
| hexadienal; 2,4-hexadienal | 637564 | value | 0.641 | 0.64 | 0.999 | update | 2026-10-23 12:16-12:27 |
| hexahydro-1h-azepine | 8119 | value | 10.4 | 1.08 | 0.104 | update | 2026-10-23 12:16-12:27 |
| hexahydrophthalic anhydride | 85689 | value | 0.0049 | 0.0009 | 0.184 | update | 2026-10-23 12:16-12:27 |
| hexanal | 6184 | manual value | 1.47 | 1.51 | 1.02 | update | 2026-10-23 12:52:34 |
| hexane; n-hexane | 8058 | value | 19.1 | 20.4 | 1.07 | update | 2026-10-23 12:16-12:27 |
| hexanedioic acid | 196 | value | 0.005 | 3.02e-08 | 6.04e-06 | update | 2026-10-23 12:16-12:27 |
| hexanoic acid | 8892 | value | 0.0197 | 0.0255 | 1.3 | update | 2026-10-23 12:16-12:27 |
| hexanol; 1-hexanol | 8103 | manual value | 1.39 | 0.127 | 0.0917 | update | 2026-10-23 12:53:10 |
| hexene; 1-hexene | 11597 | value | 22.6 | 24.5 | 1.08 | update | 2026-10-23 12:16-12:27 |
| hydrocyanic acid | 768 | value | 91.5 | 99.5 | 1.09 | update | 2026-10-23 12:16-12:27 |
| hydrogen sulfide | 402 | value | 1.89e+03 | 2.05e+03 | 1.08 | update | 2026-10-23 12:16-12:27 |
| hydroxybenzaldehyde; 2-hydroxybenzaldehyde | 6998 | value | 0.101 | 0.0791 | 0.783 | update | 2026-10-23 12:16-12:27 |
| indole | 798 | value | 0.00148 | 0.00163 | 1.1 | update | 2026-10-23 12:16-12:27 |
| iodomethane | 6328 | value | 58.8 | 54 | 0.918 | update | 2026-10-23 12:16-12:27 |
| ionol; 2,6-di-tert-butyl-4-methylphenol | 31404 | value | 0.00133 | 0.000688 | 0.519 | update | 2026-10-23 12:16-12:27 |
| isoamyl acetate; isopentyl acetate; 3-methyl-1-butanol acetate; 1-butanol, 3-methyl-, acetate | 31276 | manual value | 0.673 | 0.747 | 1.11 | update | 2026-10-23 12:55:29 |
| isoborneol | 64685 | value | 0.00485 | 0.00669 | 1.38 | update | 2026-10-23 12:16-12:27 |
| isobornyl acetate | 6950273 | value | 0.0143 |  |  | delete | 2026-10-23 12:16-12:27 |
| isobutane; 2-methylpropane | 6360 | value | 319 | 348 | 1.09 | update | 2026-10-23 12:16-12:27 |
| isobutyl acetate; acetic acid, 2-methylpropyl ester | 8038 | manual value | 2.12 | 2.37 | 1.12 | update | 2026-10-23 12:52:44 |
| isobutyl phthalate; diisobutyl phthalate | 6782 | value | 0.000301 | 6.35e-06 | 0.0211 | update | 2026-10-23 12:16-12:27 |
| isobutyric acid; 2-methylpropanoic acid; 2-methyl-propanoic acid; propanoic acid, 2-methyl- | 6590 | value | 0.204 | 0.241 | 1.18 | update | 2026-10-23 12:16-12:27 |
| isoprene; 2-methyl-1,3-butadiene; 1,3-butadiene, 2-methyl- | 6557 | value | 63.3 | 73.3 | 1.16 | update | 2026-10-23 12:16-12:27 |
| isopropenyl acetate | 7916 | value | 6.01 | 6 | 0.998 | update | 2026-10-23 12:16-12:27 |
| isopropenylbenzene | 7407 | value | 0.266 | 0.253 | 0.954 | update | 2026-10-23 12:16-12:27 |
| isopropyl acetate; acetic acid, 1-methylethyl ester | 7915 | value | 6.75 | 8.05 | 1.19 | update | 2026-10-23 12:16-12:27 |
| isopropyl formate | 12257 | value | 16.7 | 18.4 | 1.1 | update | 2026-10-23 12:16-12:27 |
| isopropylbenzene | 7406 | value | 1.2 | 0.6 | 0.5 | update | 2026-10-23 12:16-12:27 |
| jasmone; (Z)-jasmone | 1549018 | value | 0.00387 |  |  | delete | 2026-10-23 12:16-12:27 |
| l-cysteine | 5862 | value | 1.83e-07 | 2.76e-07 | 1.51 | update | 2026-10-23 12:16-12:27 |
| l-menthol | 16666 | value | 0.0707 | 0.008 | 0.113 | update | 2026-10-23 12:16-12:27 |
| l-methionine | 6137 | value | 4.01e-08 | 6.93e-08 | 1.73 | update | 2026-10-23 12:16-12:27 |
| l-proline | 145742 | value | 4.03e-09 |  |  | delete | 2026-10-23 12:16-12:27 |
| lactic acid; 2-hydroxypropanoic acid | 612 | value | 0.0108 | 0.0108 | 1.01 | update | 2026-10-23 12:16-12:27 |
| linalool | 6549 | manual value | 0.0212 | 0.021 | 0.992 | update | 2026-10-23 12:37:11 |
| linalyl acetate | 8294 | value | 0.01 | 0.0148 | 1.48 | update | 2026-10-23 12:16-12:27 |
| maleic anhydride | 7923 | value | 0.0286 | 0.0333 | 1.17 | update | 2026-10-23 12:16-12:27 |
| malonic acid | 867 | value | 0.0001 | 3.6e-07 | 0.00359 | update | 2026-10-23 12:16-12:27 |
| menthol | 1254 | value | 0.00584 | 0.00849 | 1.45 | update | 2026-10-23 12:16-12:27 |
| methanamine | 6329 | value | 306 | 353 | 1.16 | update | 2026-10-23 12:16-12:27 |
| methanethiol | 878 | value | 192 | 201 | 1.05 | update | 2026-10-23 12:16-12:27 |
| methanol; methyl alcohol | 887 | manual value | 14.2 | 16.9 | 1.19 | update | 2026-10-23 13:06:26 |
| methoxybenzene | 7519 | value | 0.471 | 0.471 | 0.999 | update | 2026-10-23 12:16-12:27 |
| methyl 2-(3-oxo-2-pent-2-enylcyclopentyl)acetate | 62388 | value | 4.49e-05 |  |  | delete | 2026-10-23 12:16-12:27 |
| methyl 2-[3-oxo-2-[(E)-pent-2-enyl]cyclopentyl]acetate | 5367719 | value | 4.49e-05 |  |  | delete | 2026-10-23 12:16-12:27 |
| methyl 2-hydroxypropanoate | 11040 | value | 0.357 | 0.467 | 1.31 | update | 2026-10-23 12:16-12:27 |
| methyl 2-methyl-2-propenoate | 6658 | value | 4.57 | 5.13 | 1.12 | update | 2026-10-23 12:16-12:27 |
| methyl 2-propenoate | 7294 | value | 8.66 | 11.5 | 1.33 | update | 2026-10-23 12:16-12:27 |
| methyl 3-oxobutanoate | 7757 | value | 0.136 | 0.119 | 0.875 | update | 2026-10-23 12:16-12:27 |
| methyl acetate; acetic acid, methyl ester | 6584 | manual value | 26.4 | 28.8 | 1.09 | update | 2026-10-23 12:49:28 |
| methyl anthranilate | 8635 | value | 0.00314 | 0.00361 | 1.15 | update | 2026-10-23 12:16-12:27 |
| methyl benzoate | 7150 | manual value | 0.0919 | 0.0507 | 0.551 | update | 2026-10-23 13:01:36 |
| methyl butyrate; methyl butanoate; butanoic acid, methyl ester | 12180 | manual value | 4.65 | 4.31 | 0.926 | update | 2026-10-23 12:50:34 |
| methyl caproate; methyl hexanoate; hexanoic acid, methyl ester | 7824 | manual value | 0.581 | 0.496 | 0.853 | update | 2026-10-23 12:55:19 |
| methyl cyanide | 6342 | value | 10.5 | 10.9 | 1.04 | update | 2026-10-23 12:16-12:27 |
| methyl decanoate | 8050 | value | 0.00447 | 0.00493 | 1.1 | update | 2026-10-23 12:16-12:27 |
| methyl formate; formic acid, methyl ester | 7865 | value | 68.4 | 78.1 | 1.14 | update | 2026-10-23 12:16-12:27 |
| methyl isothiocyanate | 11167 | value | 2.06 | 0.472 | 0.229 | update | 2026-10-23 12:16-12:27 |
| methyl jasmonate | 5281929 | value | 4.49e-05 |  |  | delete | 2026-10-23 12:16-12:27 |
| methyl pentanoate | 12206 | value | 8.4 |  |  | delete | 2026-10-23 12:16-12:27 |
| methyl phenylacetate | 7559 | value | 0.0193 | 0.0173 | 0.897 | update | 2026-10-23 12:16-12:27 |
| methyl phthalate; dimethyl phthalate | 8554 | value | 0.0229 | 0.000411 | 0.0179 | update | 2026-10-23 12:16-12:27 |
| methyl propanoate; methyl propionate | 11124 | value | 8.31 | 11.2 | 1.35 | update | 2026-10-23 12:16-12:27 |
| methyl salicylate | 4133 | manual value | 0.00719 | 0.00457 | 0.636 | update | 2026-10-23 13:01:44 |
| methyl thiocyanate | 11168 | value | 1.52 | 1.61 | 1.06 | update | 2026-10-23 12:16-12:27 |
| methylacetophenone; 1-phenyl-1-propanone | 7148 | value | 0.111 | 0.2 | 1.8 | update | 2026-10-23 12:16-12:27 |
| methylaniline; n-methylaniline | 7515 | value | 0.0778 | 0.0604 | 0.776 | update | 2026-10-23 12:16-12:27 |
| methylbutenol; 2-methyl-3-buten-2-ol | 8257 | value | 3.13 | 3.13 | 1 | update | 2026-10-23 12:16-12:27 |
| methylcyclopentane | 7296 | value | 16.7 | 18.4 | 1.1 | update | 2026-10-23 12:16-12:27 |
| methyleugenol; eugenol methyl ether | 7127 | value | 0.002 | 0.00133 | 0.667 | update | 2026-10-23 12:16-12:27 |
| methylfuran; 2-methylfuran; furan, 2-methyl- | 10797 | value | 18.9 | 18.5 | 0.979 | update | 2026-10-23 12:16-12:27 |
| methylhydrazine | 6061 | value | 5.96 | 6.64 | 1.11 | update | 2026-10-23 12:16-12:27 |
| methylnaphthalene; 1-methylnaphthalene | 7002 | value | 0.772 | 0.00893 | 0.0116 | update | 2026-10-23 12:16-12:27 |
| methylpyridine; 2-methylpyridine | 7975 | value | 1.28 | 1.49 | 1.17 | update | 2026-10-23 12:16-12:27 |
| methylpyrrole; 1-methylpyrrole | 7304 | value | 2.85 | 2 | 0.701 | update | 2026-10-23 12:16-12:27 |
| morpholine | 8083 | value | 1.04 | 1.35 | 1.3 | update | 2026-10-23 12:16-12:27 |
| myricetin | 5281672 | value | 0.0912 |  |  | delete | 2026-10-23 12:16-12:27 |
| n,n-dimethylacetamide | 31374 | value | 0.433 | 0.22 | 0.508 | update | 2026-10-23 12:16-12:27 |
| n,n-dimethylaniline | 949 | value | 0.112 | 0.0933 | 0.83 | update | 2026-10-23 12:16-12:27 |
| N,N-dimethylcyclohexanamine | 7415 | value | 0.4 | 0.4 | 1 | update | 2026-10-23 12:16-12:27 |
| n,n-dimethylformamide | 6228 | value | 0.469 | 0.503 | 1.07 | update | 2026-10-23 12:16-12:27 |
| n-butyl-n-nitroso-1-butanamine | 13542 | value | 0.00579 | 0.00625 | 1.08 | update | 2026-10-23 12:16-12:27 |
| n-butylmethylamine | 8068 | value | 6.7 | 6.73 | 1 | update | 2026-10-23 12:16-12:27 |
| n-ethylaniline | 7670 | value | 0.147 | 0.0272 | 0.185 | update | 2026-10-23 12:16-12:27 |
| n-ethylethanamine | 8021 | value | 27.8 | 31.6 | 1.14 | update | 2026-10-23 12:16-12:27 |
| n-methylformamide | 31254 | value | 0.0312 | 0.0337 | 1.08 | update | 2026-10-23 12:16-12:27 |
| n-methylmethanamine | 674 | value | 187 | 203 | 1.08 | update | 2026-10-23 12:16-12:27 |
| n-nitrosodimethylamine | 6124 | value | 1.2 | 0.36 | 0.299 | update | 2026-10-23 12:16-12:27 |
| N-propan-2-ylpropan-2-amine | 7912 | value | 11.1 | 10.6 | 0.957 | update | 2026-10-23 12:16-12:27 |
| naphthalene | 931 | value | 0.0102 | 0.0112 | 1.1 | update | 2026-10-23 12:16-12:27 |
| neoisomenthol; (+)-neoisomenthol | 19244 | value | 0.0707 |  |  | delete | 2026-10-23 12:16-12:27 |
| nicotine | 89594 | value | 0.0294 | 0.00507 | 0.172 | update | 2026-10-23 12:16-12:27 |
| nitrobenzene | 7416 | value | 0.0616 | 0.04 | 0.649 | update | 2026-10-23 12:16-12:27 |
| nitromethane | 6375 | value | 4.07 | 4.77 | 1.17 | update | 2026-10-23 12:16-12:27 |
| nonane | 8141 | value | 0.501 | 0.592 | 1.18 | update | 2026-10-23 12:16-12:27 |
| nonanone; 2-nonanone | 13187 | value | 0.148 | 0.0827 | 0.559 | update | 2026-10-23 12:16-12:27 |
| nonene; 1-nonene | 31285 | value | 0.963 | 0.72 | 0.748 | update | 2026-10-23 12:16-12:27 |
| nootkatone | 1268142 | value | 0.0004 |  |  | delete | 2026-10-23 12:16-12:27 |
| o-chloroaniline | 7240 | value | 0.0593 | 0.0272 | 0.459 | update | 2026-10-23 12:16-12:27 |
| o-cresol; 2-methylphenol | 335 | manual value | 0.0535 | 0.0319 | 0.596 | update | 2026-10-23 12:56:22 |
| octadecanoic acid | 5281 | value | 0.0888 | 9.63e-08 | 1.08e-06 | update | 2026-10-23 12:16-12:27 |
| octanal | 454 | value | 0.132 | 0.157 | 1.2 | update | 2026-10-23 12:16-12:27 |
| octane | 356 | value | 1.51 | 1.88 | 1.24 | update | 2026-10-23 12:16-12:27 |
| octanoic acid | 379 | value | 0.0448 | 0.000495 | 0.011 | update | 2026-10-23 12:16-12:27 |
| octanol; 1-octanol | 957 | manual value | 0.00954 | 0.0106 | 1.11 | update | 2026-10-23 13:00:47 |
| octene; 1-octene | 8125 | value | 2.21 | 2.32 | 1.05 | update | 2026-10-23 12:16-12:27 |
| oxalic acid | 971 | value | 0.000133 | 3.12e-05 | 0.234 | update | 2026-10-23 12:16-12:27 |
| oxirane; 1,2-epoxyethane | 6354 | value | 147 | 146 | 0.993 | update | 2026-10-23 12:16-12:27 |
| p-chloroaniline | 7812 | value | 0.00363 | 0.00947 | 2.61 | update | 2026-10-23 12:16-12:27 |
| p-cresol; 4-methylphenol | 2879 | value | 0.0147 | 0.0148 | 1 | update | 2026-10-23 12:16-12:27 |
| p-cymene; 1-isopropyl-4-methylbenzene; benzene, 1-methyl-4-(1-methylethyl)- | 7463 | value | 0.198 | 0.2 | 1.01 | update | 2026-10-23 12:16-12:27 |
| pentadecane | 12391 | value | 0.000557 | 0.000656 | 1.18 | update | 2026-10-23 12:16-12:27 |
| pentadiene; 1,3-pentadiene | 62204 | value | 54.2 | 54 | 0.997 | update | 2026-10-23 12:16-12:27 |
| pentanal; valeraldehyde | 8063 | manual value | 3.99 | 3.43 | 0.861 | update | 2026-10-23 12:50:27 |
| pentane | 8003 | value | 59.9 | 68.5 | 1.14 | update | 2026-10-23 12:16-12:27 |
| pentanedioic acid | 743 | value | 0.0667 | 3.84e-07 | 5.76e-06 | update | 2026-10-23 12:16-12:27 |
| pentanenitrile | 8061 | value | 1.15 | 0.973 | 0.844 | update | 2026-10-23 12:16-12:27 |
| pentanol; 1-pentanol | 6276 | manual value | 0.39 | 0.293 | 0.752 | update | 2026-10-23 12:51:01 |
| pentene; 1-pentene | 8004 | value | 46.4 | 59.1 | 1.28 | update | 2026-10-23 12:16-12:27 |
| pentyl acetate | 12348 | manual value | 0.583 | 0.558 | 0.958 | update | 2026-10-23 12:55:20 |
| phenanthrene | 995 | value | 0.0445 | 1.61e-05 | 0.000363 | update | 2026-10-23 12:16-12:27 |
| phenethyl propanoate; phenethyl propionate | 31225 | value | 0.00673 | 0.00667 | 0.99 | update | 2026-10-23 12:16-12:27 |
| phenethylamine; phenylethylamine; benzeneethanamine | 1001 | value | 0.0348 | 0.031 | 0.89 | update | 2026-10-23 12:16-12:27 |
| phenol | 996 | value | 0.0492 | 0.0467 | 0.95 | update | 2026-10-23 12:16-12:27 |
| phenylacetonitrile | 8794 | value | 0.0113 | 0.0119 | 1.05 | update | 2026-10-23 12:16-12:27 |
| phthalic anhydride | 6811 | value | 0.000113 | 6.89e-05 | 0.611 | update | 2026-10-23 12:16-12:27 |
| pinane | 10129 | value | 0.416 | 0.688 | 1.66 | update | 2026-10-23 12:16-12:27 |
| piperidine | 8082 | value | 4.91 | 4.28 | 0.873 | update | 2026-10-23 12:16-12:27 |
| piperonal | 8438 | value | 0.0673 | 0.00133 | 0.0198 | update | 2026-10-23 12:16-12:27 |
| propanal | 527 | manual value | 35.7 | 42.3 | 1.18 | update | 2026-10-23 12:49:26 |
| propane | 6334 | value | 961 | 953 | 0.992 | update | 2026-10-23 12:16-12:27 |
| propanenitrile | 7854 | value | 5.53 | 6.32 | 1.14 | update | 2026-10-23 12:16-12:27 |
| propanoic acid; propionic acid | 1032 | value | 0.703 | 0.471 | 0.669 | update | 2026-10-23 12:16-12:27 |
| propyl acetate; n-propyl acetate; acetic acid, propyl ester | 7997 | manual value | 4.75 | 4.79 | 1.01 | update | 2026-10-23 12:50:34 |
| propyl glycol; 2-propoxyethanol; 2-propyloxyethanol; ethylene glycol monopropyl ether | 17756 | value | 0.311 | 0.258 | 0.831 | update | 2026-10-23 12:16-12:27 |
| propylbenzene | 7668 | value | 0.638 | 0.456 | 0.714 | update | 2026-10-23 12:16-12:27 |
| propylene carbonate | 7924 | value | 0.00567 | 0.006 | 1.06 | update | 2026-10-23 12:16-12:27 |
| psoralen | 6199 | value | 7.2e-07 |  |  | delete | 2026-10-23 12:16-12:27 |
| putrescine; 1,4-butanediamine | 1045 | value | 0.43 | 0.549 | 1.28 | update | 2026-10-23 12:16-12:27 |
| pyrene | 31423 | value | 0.116 | 6e-07 | 5.19e-06 | update | 2026-10-23 12:16-12:27 |
| pyridine | 1049 | value | 2.41 | 2.77 | 1.15 | update | 2026-10-23 12:16-12:27 |
| pyrrole | 8027 | value | 1.11 | 1.11 | 0.999 | update | 2026-10-23 12:16-12:27 |
| pyrrolidine | 31268 | value | 6.17 | 8.36 | 1.35 | update | 2026-10-23 12:16-12:27 |
| quercetin | 5280343 | value | 0.0375 |  |  | delete | 2026-10-23 12:16-12:27 |
| quinoline | 7047 | value | 0.0393 | 0.008 | 0.203 | update | 2026-10-23 12:16-12:27 |
| sec-amyl acetate; 2-pentanol acetate; 2-pentanol, acetate; 1-methylbutyl acetate | 12278 | value | 0.994 | 0.932 | 0.937 | update | 2026-10-23 12:16-12:27 |
| sec-butyl acetate | 7758 | value | 2.57 | 1.95 | 0.76 | update | 2026-10-23 12:16-12:27 |
| squalene | 638072 | value | 8.4e-07 |  |  | delete | 2026-10-23 12:16-12:27 |
| stigmasterol | 5280794 | value | 0.0216 |  |  | delete | 2026-10-23 12:16-12:27 |
| styrene; vinylbenzene | 7501 | value | 0.829 | 0.853 | 1.03 | update | 2026-10-23 12:16-12:27 |
| sulfur dioxide | 1119 | value | 341 | 327 | 0.961 | update | 2026-10-23 12:16-12:27 |
| taurine | 1123 | value | 1.15e-05 | 2.29e-05 | 2 | update | 2026-10-23 12:16-12:27 |
| terpineol; alpha-terpineol; .alpha.-terpineol | 17100 | manual value | 0.00549 | 0.00564 | 1.03 | update | 2026-10-23 12:37:25 |
| tert-butyl ethyl ether | 12512 | value | 8 | 12.8 | 1.6 | update | 2026-10-23 12:16-12:27 |
| tert-butyl methyl ether | 15413 | value | 31 | 32.7 | 1.05 | update | 2026-10-23 12:16-12:27 |
| tetrabromomethane | 11205 | value | 3.59 | 0.096 | 0.0268 | update | 2026-10-23 12:16-12:27 |
| tetrachloroethene | 31373 | value | 2.11 | 2.47 | 1.17 | update | 2026-10-23 12:16-12:27 |
| tetrachloromethane | 5943 | value | 13.7 | 15.3 | 1.12 | update | 2026-10-23 12:16-12:27 |
| tetradecane | 12389 | value | 0.0456 | 0.00177 | 0.0389 | update | 2026-10-23 12:16-12:27 |
| tetradecanoic acid | 11005 | value | 0.0444 | 1.87e-07 | 4.2e-06 | update | 2026-10-23 12:16-12:27 |
| tetradecene; 1-tetradecene | 14260 | value | 0.153 | 0.002 | 0.0131 | update | 2026-10-23 12:16-12:27 |
| tetrahydrofuran | 8028 | value | 19.5 | 21.6 | 1.11 | update | 2026-10-23 12:16-12:27 |
| tetrahydrofurfuryl alcohol | 7360 | value | 0.133 | 0.146 | 1.1 | update | 2026-10-23 12:16-12:27 |
| tetrahydrothiophene | 1127 | value | 2.43 | 2.4 | 0.989 | update | 2026-10-23 12:16-12:27 |
| tetramethylurea | 12437 | value | 1.85 |  |  | delete | 2026-10-23 12:16-12:27 |
| thiophene | 8030 | value | 8.64 | 10.6 | 1.23 | update | 2026-10-23 12:16-12:27 |
| thiourea | 2723790 | value | 0.000187 | 1.88e-08 | 0.000101 | update | 2026-10-23 12:16-12:27 |
| thymol | 6989 | value | 0.0371 | 0.00213 | 0.0575 | update | 2026-10-23 12:16-12:27 |
| toluene; methylbenzene; benzene, methyl- | 1140 | value | 3.27 | 3.79 | 1.16 | update | 2026-10-23 12:16-12:27 |
| trans-(2R,5S)-5-methyl-2-propan-2-ylcyclohexan-1-one | 443159 | value | 0.036 |  |  | delete | 2026-10-23 12:16-12:27 |
| trans-resveratrol | 445154 | value | 1.65e-10 |  |  | delete | 2026-10-23 12:16-12:27 |
| triacetin; 1,2,3-propanetriol triacetate | 5541 | value | 0.0336 | 0.00033 | 0.00984 | update | 2026-10-23 12:16-12:27 |
| tribromomethane | 5558 | value | 0.694 | 0.67 | 0.965 | update | 2026-10-23 12:16-12:27 |
| tributyl acetylcitrate | 6505 | value | 0.0445 | 6.07e-07 | 1.36e-05 | update | 2026-10-23 12:16-12:27 |
| tributylamine | 7622 | value | 0.897 | 0.0125 | 0.0139 | update | 2026-10-23 12:16-12:27 |
| trichloroacetic acid | 6421 | value | 0.221 | 0.008 | 0.0362 | update | 2026-10-23 12:16-12:27 |
| trichlorobenzene | 6895 | value | 0.0517 | 0.0244 | 0.473 | update | 2026-10-23 12:16-12:27 |
| trichloroethane; 1,1,1-trichloroethane | 6278 | value | 14.9 | 16.5 | 1.11 | update | 2026-10-23 12:16-12:27 |
| trichloromethane | 6212 | value | 47.7 | 26.3 | 0.551 | update | 2026-10-23 12:16-12:27 |
| tridecane | 12388 | value | 0.00583 | 0.00623 | 1.07 | update | 2026-10-23 12:16-12:27 |
| tridecanone; 2-tridecanone | 11622 | value | 0.133 |  |  | delete | 2026-10-23 12:16-12:27 |
| triethyl citrate | 6506 | value | 0.0446 | 0.0003 | 0.00672 | update | 2026-10-23 12:16-12:27 |
| triethylamine | 8471 | value | 7.32 | 7.61 | 1.04 | update | 2026-10-23 12:16-12:27 |
| trimethylamine | 1146 | value | 198 | 215 | 1.09 | update | 2026-10-23 12:16-12:27 |
| trimethylbenzene; 1,2,3-trimethylbenzene | 10686 | value | 0.191 | 0.225 | 1.18 | update | 2026-10-23 12:16-12:27 |
| tryptophan | 6305 | value | 2.8e-10 |  |  | delete | 2026-10-23 12:16-12:27 |
| undecane | 14257 | value | 0.0548 | 0.0549 | 1 | update | 2026-10-23 12:16-12:27 |
| undecene; 1-undecene | 13190 | value | 0.0655 | 0.0657 | 1 | update | 2026-10-23 12:16-12:27 |
| ursolic acid | 64945 | value | 0.0465 |  |  | delete | 2026-10-23 12:16-12:27 |
| valeric acid; pentanoic acid | 7991 | value | 0.0512 | 0.0261 | 0.51 | update | 2026-10-23 12:16-12:27 |
| vanillin | 1183 | value | 0.0445 | 1.57e-05 | 0.000354 | update | 2026-10-23 12:16-12:27 |
| vinyl acetate; acetic acid ethenyl ester | 7904 | value | 12.3 | 11.7 | 0.952 | update | 2026-10-23 12:16-12:27 |
| vinylcyclohexene | 17497 | null |  | 2 |  | update | 2026-10-23 12:57:32 |
| vinylpyridine | 7521 | value | 1.33 |  |  | delete | 2026-10-23 12:16-12:27 |
