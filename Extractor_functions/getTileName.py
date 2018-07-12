from landsat_index import ConvertToWRS_Landsat
from sentinel_index import ConvertToWRS_Sentinel
from os.path import join
import os, sys
import pprint
import numpy as np
import subprocess

sys.path.append(join(os.path.dirname(os.path.realpath(__file__)),".."))
printer = pprint.PrettyPrinter(indent=3)

'''
For Landsat Data
'''
def get_landsat_tile(points):
    tile_dict = {}
    to_WRS = ConvertToWRS_Landsat()
    for (lat,lon) in points:
        result = to_WRS.get_wrs(lat,lon)
        for i in range(0,len(result)):
            path = '0' + str(result[i]['path'])
            row = '0' + str(result[i]['row'])
            key = path + '-' + row
            tile_dict.setdefault(key,[]).append((lat,lon))
    return tile_dict

'''
For Sentinel Data
'''
def get_sentinel_tile(points):
    tile_dict = {}
    to_WRS = ConvertToWRS_Sentinel()
    for (lat,lon) in points:
        result = to_WRS.get_wrs(lat,lon)
        for key in result:
            tile_dict.setdefault(key,[]).append((lat,lon))
    return tile_dict

'''
For MODIS Data
'''
def get_modis_tile(points):
    tile_dict = {}
    for pt in points:
        try:
            tileName = subprocess.check_output(['extractor_functions/tilemap3', str(pt[0]), str(pt[1])])
            tileName = tileName.decode("utf-8")
            tileSub = tileName.split('#')[1]
            tileName = 'h' + tileSub[0: -3].zfill(2) + 'v' + tileSub[-3:-1]
        except Exception as e:
            print("Call MODIS2TILE program error: " + str(pt) + "\n")
            continue
        if tileName is None:
            continue
        if tileName in tile_dict.keys():
            tile_dict[tileName].append(pt)
        else:
            tile_dict[tileName] = [pt]
    return tile_dict

if __name__ == '__main__':
    print('begin!')
    points = ((np.load('2016TT_points.npz')['arr_0']).tolist())['Corn']
    printer.pprint(get_landsat_tile(points).keys())
    # printer.pprint(get_sentinel_tile(points))
    # printer.pprint(get_modis_tile(points))
    print('finish!')
