# use conda environment with package pip install "cdsapi>=0.7.7" ('test' for Sarah Grandke)
# need file '.cdsapirc' in home folder with:
#   url: https://ads.atmosphere.copernicus.eu/api
#   key: cad840ed-3c1f-4249-9901-b2f93d9886b5

# you can unzip the resulting files all at once using the following command while in the folder with the files.
# for z in cams_co2_*; do unzip $z; done

import cdsapi

for l in ["2020","2021","2022"]:
    if l=="2020": list= ["10", "11", "12"]
    if l=="2021": list= ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    if l=="2022": list= ["01", "02", "03"]
    for k in list:
        dataset = "cams-global-greenhouse-gas-inversion"
        request = {
            "variable": "carbon_dioxide",
            "quantity": "concentration",
            "input_observations": "satellite",
            "time_aggregation": "instantaneous",
            "version": "latest",
            "year": [l],
            "month": [k]
        }

        client = cdsapi.Client()
        client.retrieve(dataset, request).download(f"/work/bb1170/RUN/b383736/data/CAMS/2020-2022/satellite/cams_co2_{l}_{k}")

