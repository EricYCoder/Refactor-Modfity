"""
function to extract LST data (MODIS) given latitude, longitude and a time period
"""

import os, sys, time, pprint
from osgeo import gdal
from osgeo import osr
from os.path import join
import numpy as np
from datetime import date
from osgeo.gdalconst import *
from datetime import datetime
import subprocess
import settings

sys.path.append(join(os.path.dirname(os.path.realpath(__file__)),".."))
printer = pprint.PrettyPrinter(indent=3)

'''
function to get hdf files based on tileName and Date range
'''
def create_tar_hdf(data_path, tile_Key, startDate, endDate):
    try:
        startDate = datetime.strptime(startDate, "%Y%m%d")
        endDate = datetime.strptime(endDate, "%Y%m%d")
        if startDate > endDate:
            print("startDate > endDate")
            return[]

        hdf_list = []
        for MOD11Folder in os.listdir(data_path):
            MOD11Date = datetime.strptime(MOD11Folder, "%Y.%m.%d")
            if startDate <= MOD11Date <= endDate:
                for MOD11_file in os.listdir(join(data_path, MOD11Folder)):
                    if (MOD11_file.split(".")[-1] == "hdf") and (MOD11_file.split(".")[2] in [tile_Key]):
                        hdf_list.append(join(data_path, MOD11Folder, MOD11_file))
            else:
                continue

        return(hdf_list)

    except:
        return[]

'''
function to convert lat and lng to projection coordinate
'''
def latlng2geo(dataset, lat, lng):
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()

    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lng, lat)
    return coords[:2]

'''
function to convert projection coordinate to image coordinate
'''
def geo2imagexy(dataset, proj_x, proj_y):
    trans = dataset.GetGeoTransform()
    dTemp = trans[1] * trans[5] - trans[2] * trans[4]
    img_x = (int)((trans[5] * (proj_x - trans[0]) - trans[2] * (proj_y - trans[3])) / dTemp)
    img_y = (int)((trans[1] * (proj_y - trans[3]) - trans[4] * (proj_x - trans[0])) / dTemp)
    return (img_x, img_y)

'''
function to convert date(YYYYMMDD) to DOY
'''
def date2DOY(file_date):
    try:
        file_date = datetime.strptime(file_date, "%Y%m%d")
        days = str(file_date.year).zfill(4) + str(file_date.timetuple().tm_yday).zfill(3)
        return days
    except:
        return None

'''
Get tilename based on lat and lng
'''
def get_tile_name(points):
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

        if tileName == None:
            continue
        if tileName in tile_dict.keys():
            tile_dict[tileName].append(pt)
        else:
            tile_dict[tileName] = [pt]
    return tile_dict

'''
Get band value based on image coordinate
'''
def get_band_value(bandraster, px, py, cloud=False):
    result = bandraster[py][px]
    if cloud is True:
        value = round(result, 4)
    elif 7500 <= result <= 65535:
        value = round(result*0.02, 4)
    else:
        value = -1
    return value

'''
Get image coordinate based on lat/lon
'''
def point_boundary_is_valid(bandfile, lat, lon):
    xsize = bandfile.RasterXSize
    ysize = bandfile.RasterYSize
    x, y = latlng2geo(bandfile, lat, lon) # (lat, lon) to (x, y)
    px, py = geo2imagexy(bandfile, x, y)  # (x, y) to (row, col)
    px = int(px)
    py = int(py)
    if px > xsize or px < 0 or py > ysize or py < 0:
        printer.pprint("the point is out of the image range:" + str(lat) + ', ' + str(lon))
        return False, px, py
    else:
        return True, px, py

'''
function to extract MODIS LST data given latitude, longitude and a period of time
input of the function:
points: coordinate (lat, lon); float, a list of tuple
startDate: start time, should be in the format of 'YYYYMMDD' [string]
endDate: end time, should be in the format of "YYYYMMDD" [string]
output of the function:
a list of data with format[{coordinate: value, 'MOD11A1': { day/night: [(time1,value1),(time2,value2),....]}, 'MYD11A1': { day/night: [(time1,value1),(time2,value2),....]}},...]
---------------Valid LST Data-----------------
7500 <= LST <= 65535 (scale is 0.02)
QC == 0 or QC & 0x000F == 1
--------------Invalid LST Data----------------
LST = -1
Others(eg. QC == 2/3)
'''

def extract_MODIS_LST(points, startDate, endDate):
    #Remove duplicate points
    #points = list(set(points))

    #Results Initialization
    MODIS_LST_res = []
    for pt in points:
        MODIS_LST_record = {}
        #Add coordinate key
        MODIS_LST_record['coordinate'] = (pt[0], pt[1])
        MODIS_LST_record['MOD11A1'] = []
        MODIS_LST_record['MYD11A1'] = []
        MODIS_LST_res.append(MODIS_LST_record)

    #Get tile name based on lat/lon
    time00 = time.time()
    tile_dict = get_tile_name(points)
    if len(tile_dict) == 0:
        raise Exception('There is no consitent tile')
        return
    else:
        pass
    print('Get tile time:' + str(time.time() - time00))

    #MODIS LST folders
    MOD_data_path = join(os.path.expanduser('~'), settings.MOD11A1_PATH)
    MYD_data_path = join(os.path.expanduser('~'), settings.MYD11A1_PATH)

    #Get Data
    for tile_Key in tile_dict.keys():
        hdf_files = []
        hdf_files = hdf_files + create_tar_hdf(MOD_data_path, tile_Key, startDate, endDate)
        hdf_files = hdf_files + create_tar_hdf(MYD_data_path, tile_Key, startDate, endDate)

        for file_path in hdf_files:
            try:
                hdf_ds = gdal.Open(file_path, gdal.GA_ReadOnly)
            except:
                print('Unable to open file: \n' + file_path + '\nSkipping\n')
                continue

            #Compare file info and folder info
            path_info = file_path.split('/')
            file_date = path_info[-2].replace('.', '')
            file_type = path_info[-3].split('.')[0]
            file_name = os.path.basename(file_path)
            file_info = file_name.split('.')
            if file_type != file_info[0] or date2DOY(file_date) != file_info[1][1:]:
                print('File info unmatch: \n' + file_path + '\nSkipping\n')
                continue

            #read dataset
            #Get Day LST and QC
            #Get Day QC Dataset
            day_lst_ds = gdal.Open(hdf_ds.GetSubDatasets()[0][0])
            day_lst_rst = day_lst_ds.GetRasterBand(1).ReadAsArray()
            day_qc_ds = gdal.Open(hdf_ds.GetSubDatasets()[1][0])
            day_qc_rst = day_qc_ds.GetRasterBand(1).ReadAsArray()

            #Get Night LST and QC
            #Get Night QC Value
            night_lst_ds = gdal.Open(hdf_ds.GetSubDatasets()[4][0])
            night_lst_rst = night_lst_ds.GetRasterBand(1).ReadAsArray()
            night_qc_ds = gdal.Open(hdf_ds.GetSubDatasets()[5][0])
            night_qc_rst = night_qc_ds.GetRasterBand(1).ReadAsArray()

            for location in tile_dict[tile_Key]:
                latitude = location[0]
                longitude = location[1]

                #Get img coords(from proj coord) and check if data is valid
                is_valid, px, py = point_boundary_is_valid(day_qc_ds, latitude, longitude)
                if not is_valid:
                    continue
                
                #Get Day QC Value
                day_lst_value = -1
                day_qc_value = get_band_value(day_qc_rst, px, py, True)
                if day_qc_value == 0 or day_qc_value & 0x000F == 1:
                    #Get Day LST Value
                    day_lst_value = get_band_value(day_lst_rst, px, py)

                #Get Night QC Value
                night_lst_value = -1
                night_qc_value = get_band_value(night_qc_rst, px, py, True)
                if night_qc_value == 0 or night_qc_value & 0x000F == 1:
                    #Get Night LST Value
                    night_lst_value = get_band_value(night_lst_rst, px, py)

                #Add to results
                point_exist = False

                #if the point exist
                for MODIS_LST_record in MODIS_LST_res:
                    if list(MODIS_LST_record['coordinate']) == list(location):
                        point_exist = True
                        try:
                            MODIS_LST_record[file_type].append((file_date, float('%.2f'%day_lst_value), float('%.2f'%night_lst_value)))
                            print("Extract finish in :" + file_path)
                            break
                        except KeyError:
                            pass

                if point_exist == False:
                    #if the point does not exist
                    print("The point has not been initialized: " + str(location) + "\n")
                    MODIS_LST_record = {}
                    #Add coordinate key
                    MODIS_LST_record['coordinate'] = (latitude, longitude)
                    MODIS_LST_record['MOD11A1'] = []
                    MODIS_LST_record['MYD11A1'] = []
                    #Add LST values
                    MODIS_LST_record[file_type] = [(file_date, float('%.2f'%day_lst_value), float('%.2f'%night_lst_value))]
                    MODIS_LST_res.append(MODIS_LST_record)


    #sort the results by Date
    for MODIS_LST_record in MODIS_LST_res:
        if 'MOD11A1' in MODIS_LST_record.keys():
            MODIS_LST_record['MOD11A1'].sort()
        if 'MYD11A1' in MODIS_LST_record.keys():
            MODIS_LST_record['MYD11A1'].sort()

    return (MODIS_LST_res)


# test main function
if __name__ == '__main__':
    time0 = time.time()
    # points = [(47.729849, -119.540372), (43.79562, -111.466202), (42.034932, -107.947733), (30.732852, -104.768943), (42.969349, -109.885573)]
    points = list(np.load('lonlatArray.npz')['arr_0'])
    startDate = "20100401"
    endDate = "20101001"
    MODIS_LST_res = extract_MODIS_LST(points, startDate, endDate)
    printer.pprint(MODIS_LST_res)
    print('Overall Time:' + str(time.time() - time0))
    print("over")
