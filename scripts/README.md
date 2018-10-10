# Aldy Database YAMLs: Preparation scripts

The purpose of `get_genes.py` is to process CYP databases in old Karolinska HTML format
and to generate Aldy-compatible YAML files.

### Requirements

- Python 2.7 (Python 3 is not supported)
- The following Python libraries:
       
  ```
  pip install beautifulsoup4 myvariant biopython lxml logbook
  ```

- BLAT compiled with g++. clang-compiled BLAT will not work! 
  
  > macOS only: a proper BLAT executable is available in the working directory.

- Working internet connection for NCBI and myVariant queries.

### How to run?

Just do

```
python2 get_genes.py
```

and wait. YML files will be generated in `output` directory.


### FAQ

1. _This script only generates YML for 6 genes. What about other 4?_
   
   Two genes (_CYP4F2_ and _CYP2C8_) are not supported anymore due to the 
   NCBI database changes. _DPYD_ and _TMPT_ databases were constructed manually.
   
2. _I have found bugs and problems with the scripts--- what should I do?_
   
   I would suggest manual YML editing. I haven't used these scripts in a while, and plan to retire them soon. 