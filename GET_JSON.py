    # -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 21:59:56 2025

@author: itai
"""

import requests
import json
import configparser

# Space-Track API credentials
#---------------------------------------------read user acces ----------------------------------------------
# Use configparser package to pull in the ini file (pip install configparser)
config = configparser.ConfigParser()
config.read("./UserInit.ini")
username = config.get("configuration","username")
password = config.get("configuration","password")
#configOut = config.get("configuration","output")
siteCred = {'identity': username, 'password': password}


# Base URL and TLE endpoint
uriBase = "https://www.space-track.org"
requestLogin = "/ajaxauth/login"
requestTLELatest = "/basicspacedata/query/class/tle_latest/ORDINAL/1/format/json"  # Query for latest TLEs

# Log in to Space-Track.org
with requests.Session() as session:
    # Authenticate
    credentials = {'identity': username, 'password': password}
    resp = session.post(uriBase + requestLogin, data=credentials)
    
    if resp.status_code == 200:
        print("Login successful.")
    else:
        print("Login failed:", resp.text)
        exit()

    # Fetch latest TLE data
    resp = session.get(uriBase + requestTLELatest)
    
    if resp.status_code == 200:
        print("TLE data retrieved successfully.")
        # Save the JSON response to 'celes_tle.json'
        with open("celes_tle.json", "w") as file:
            json.dump(resp.json(), file, indent=4)
        print("TLE data saved to 'celes_tle.json'.")
    else:
        print("Failed to retrieve TLE data:", resp.text)




# # Base URL and endpoints
# uriBase = "https://www.space-track.org"
# requestLogin = "/ajaxauth/login"
# requestSatcat = "/basicspacedata/query/class/satcat/format/json"  # Query to fetch all satellites

# # Log in to Space-Track.org
# with requests.Session() as session:
#     # Authenticate
#     credentials = {'identity': username, 'password': password}
#     resp = session.post(uriBase + requestLogin, data=credentials)
    
#     if resp.status_code == 200:
#         print("Login successful.")
#     else:
#         print("Login failed:", resp.text)
#         exit()

#     # Fetch satellite data
#     resp = session.get(uriBase + requestSatcat)
    
#     if resp.status_code == 200:
#         print("Satellite data retrieved successfully.")
#         # Save the JSON response to a file
#         with open("celes_tle.json", "w") as file:
#             json.dump(resp.json(), file, indent=4)
#         print("Satellite data saved to 'celes_tle.json'.")
#     else:
#         print("Failed to retrieve satellite data:", resp.text)
