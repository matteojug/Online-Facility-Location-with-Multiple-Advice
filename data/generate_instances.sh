pip3 install pyproj

python3 create_synth.py

# Add country informations to checkins
python3 enrich_loc.py raw/loc-gowalla_totalCheckins.txt
python3 enrich_loc.py raw/loc-brightkite_totalCheckins.txt
mkdir final_days
python3 create_checkins.py

mkdir final_uber
python3 create_uber.py