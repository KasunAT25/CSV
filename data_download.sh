# chmod +x data_download.sh
# ./data_download.sh

echo "downloading data ..."
cd data
for d in covid genome osm fb
do
wget --no-check-certificate -nc https://www.cse.cuhk.edu.hk/mlsys/gre/$d
done
