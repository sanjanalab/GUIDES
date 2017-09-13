pip install -r requirements.txt
deactivate
npm install
cd static
rm -rf data
wget -r --no-parent https://www.dropbox.com/s/sp59h2m6vjrsr1c/GUIDES_DATA.zip?dl=0
cd ..
