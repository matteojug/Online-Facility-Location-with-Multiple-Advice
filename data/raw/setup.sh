wget https://snap.stanford.edu/data/loc-brightkite_totalCheckins.txt.gz
wget https://snap.stanford.edu/data/loc-gowalla_totalCheckins.txt.gz
wget http://download.geonames.org/export/dump/cities500.zip

# Download uber from https://www.kaggle.com/fivethirtyeight/uber-pickups-in-new-york-city

gzip -d loc-brightkite_totalCheckins.txt.gz
gzip -d loc-gowalla_totalCheckins.txt.gz
unzip cities500.zip
