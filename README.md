# FinderChart

Tool to make finder charts using ZTF science images or reference images.

Usage: 
get_finder(RA, Dec, SourceName)

Generates a finder chart image, and prints starlist for DBSP or Keck

Borrows heavily from Nadia's PS1 finder chart code,
https://github.com/nblago/utils/blob/master/finder_chart.py

and makes use of ztfquery,
https://github.com/MickaelRigault/ztfquery
