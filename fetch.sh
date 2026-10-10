#!/bin/bash
# Download every primary source into rawdata/, unmodified.  The executable record of
# provenance; run from an empty checkout before exportweb.py / validate.py.
set -e
cd "$(dirname "$0")"
mkdir -p rawdata/draine rawdata/dl01 rawdata/ld01 rawdata/mmp83
# Draine (2003) dielectric functions and Draine's Q tables (used only for checks).
# His server refuses requests without a curl-like user agent.
for f in callindex.out_silD03 callindex.out_CpaD03_0.01 callindex.out_CpaD03_0.10 \
         callindex.out_CpeD03_0.01 callindex.out_CpeD03_0.10 Gra_81.gz suvSil_81.gz; do
  curl -sfL -A "curl/8.7.1" -o rawdata/draine/$f https://www.astro.princeton.edu/~draine/dust/diel/$f
done
# Draine & Li (2001, ApJ 551, 807) and Li & Draine (2001, ApJ 554, 778): arXiv sources
curl -sfL -o rawdata/dl01/src.tar.gz https://arxiv.org/e-print/astro-ph/0011318
curl -sfL -o rawdata/ld01/src.tar.gz https://arxiv.org/e-print/astro-ph/0011319
(cd rawdata/dl01 && tar xzf src.tar.gz)
(cd rawdata/ld01 && tar xzf src.tar.gz)
# Mathis, Mezger & Panagia (1983, A&A 128, 212): ADS scan (image only; Table A3 is
# transcribed in data/mmp83_tableA3.txt).  The ADS article server times out now and then.
for i in 1 2 3 4 5; do
  curl -sfL --max-time 90 -A "Mozilla/5.0" -o rawdata/mmp83/mmp83.pdf \
    "https://articles.adsabs.harvard.edu/pdf/1983A&A...128..212M" && file rawdata/mmp83/mmp83.pdf | grep -q PDF && break
  sleep 5
done
file rawdata/mmp83/mmp83.pdf | grep -q PDF || { echo "MMP83 scan not fetched"; exit 1; }
# Cardelli, Clayton & Mathis (1989, ApJ 345, 245): ADS scan (Eqs. 2-5 and Table 3)
mkdir -p rawdata/ccm89
for i in 1 2 3 4 5; do
  curl -sfL --max-time 90 -A "Mozilla/5.0" -o rawdata/ccm89/ccm89.pdf \
    "https://articles.adsabs.harvard.edu/pdf/1989ApJ...345..245C" && file rawdata/ccm89/ccm89.pdf | grep -q PDF && break
  sleep 5
done
echo "rawdata complete"
