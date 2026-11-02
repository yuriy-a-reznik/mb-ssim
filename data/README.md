# data/

Dataset index files go here. They are not committed because they embed the
absolute path of your copy of each database. Build them with

    python3 scripts/make_index.py cidiq /path/to/CID-IQ --distance=50  -o data/cidiq_50cm.csv
    python3 scripts/make_index.py cidiq /path/to/CID-IQ --distance=100 -o data/cidiq_100cm.csv
    python3 scripts/make_index.py live  /path/to/databaserelease2       -o data/live_r2_realigned.csv

See docs/DATASETS.md for where to obtain the databases and the viewing
geometry behind each default `--ppd`.
