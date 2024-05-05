import urllib.request as req
import os.path
import json
import sys
API_ROOT = "https://api.synopticdata.com/v2/"
API_TOKEN = "f2f69390151c4a4da5d51cc756ffa6f8"
# let's get some latest data
api_request_url = os.path.join(API_ROOT, "stations/latest")
# the built-in library requires us to specify the URL parameters directly
api_request_url += "?bbox=76,8,81,14"
api_request_url += "&vars=air_temp"
api_request_url += "&token="+API_TOKEN
print(api_request_url)
# Note, .format is a python method available to put values into strings
# then we make the request
response = req.urlopen(api_request_url)
api_text_data = response.read()
print(api_text_data)
