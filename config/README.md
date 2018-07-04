## Config

This folder contains the configuration files necessary for running the demographic models M1-M6 mentioned in the [bioRxiv](http://dx.doi.org/10.1101/077958) draft and for converting files in `pseudo_MS` format into the `cbSFS` (also used in the draft). 

### Running models M1-M6
Here we are assuming 2 diploid samples/population (which is also applicable to 4 haploid samples/population), block sizes of 2Kb and the working directory is that in which the ABLE binary is situated. The corresponding command lines for launching these models have been provided below.

##### N.B. In order to run these very same command lines on the 500bp data provided in the `data` folder, you need to replace `-r tbi2 2001` by `-r tbi2 501` and change the filename for the `datafile` keyword in the respective config file.


#### M1
```
./ABLE 8 xxx -t tbi1 -r tbi2 2001 -I 2 4 4 -ej tbi3 1 2 -T config_M1.txt
```

#### M2
```
./ABLE 8 xxx -t tbi1 -r tbi2 2001 -I 2 4 4 -n 1 tbi4 -n 2 tbi5 -ej tbi3 1 2 -T config_M2.txt
```

#### M3
```
./ABLE 8 xxx -t tbi1 -r tbi2 2001 -I 2 4 4 -n 1 tbi4 -n 2 tbi5 -g 1 tbi6 -g 2 tbi7 -ej tbi3 1 2 -T config_M3.txt
```

#### M4
```
./ABLE 8 xxx -t tbi1 -r tbi2 2001 -I 2 4 4 -n 1 tbi4 -n 2 tbi5 -m 1 2 tbi6 -m 2 1 tbi7 -ej tbi3 1 2 -T config_M4.txt
```

#### M5
```
./ABLE 8 xxx -t tbi1 -r tbi2 2001 -I 2 4 4 -n 1 tbi4 -n 2 tbi5 -em tbi8 1 2 tbi6 -em tbi8 2 1 tbi7 -ej tbi3 1 2 -T config_M5.txt
```

#### M6
```
./ABLE 8 xxx -t tbi1 -r tbi2 2001 -I 2 4 4 -n 1 tbi4 -n 2 tbi5 -es tbi8 1 tbi6 -ej tbi8 3 2 -es tbi8 2 tbi7 -ej tbi8 4 1 -ej tbi3 1 2 -T config_M6.txt
```

### Converting files from `pseudo_MS` format to the `bSFS`
The following command lines along with the provided config files can be used to produce the `bSFS` input used in the paper (which is also available in the `data` folder).

```
./ABLE config_bSFS_Pongo_500bp.txt
./ABLE config_bSFS_Pongo_2kb.txt
```

### Converting files from `pseudo_MS` format to the `cbSFS` (composite `bSFS`)
The following command lines along with the provided config files can be used to produce the `cbSFS` input used in the paper (which is also available in the `data` folder).

```
./ABLE config_cbSFS_Pongo_500bp.txt
./ABLE config_cbSFS_Pongo_2kb.txt
```

