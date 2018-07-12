"""
function to extract LST data (MODIS) given latitude, longitude and a time period
"""

import os
import sys
import pprint
from osgeo import gdal
from os.path import join
import numpy as np
from datetime import datetime
import settings
import geo_functions

sys.path.append(join(os.path.dirname(os.path.realpath(__file__)), ".."))
printer = pprint.PrettyPrinter(indent=3)

"""
function to get hdf files based on tileName and Date range
"""


def create_tar_hdf(data_path, tile_Key, startDate, endDate):
    startDate = datetime.strptime(startDate, "%Y%m%d")
    endDate = datetime.strptime(endDate, "%Y%m%d")
    if startDate > endDate:
        print("startDate > endDate")
        return []

    hdf_list = []
    for MOD11Folder in os.listdir(data_path):
        MOD11Date = datetime.strptime(MOD11Folder, "%Y.%m.%d")
        if startDate <= MOD11Date <= endDate:
            for MOD11_file in os.listdir(join(data_path, MOD11Folder)):
                if (MOD11_file.split(".")[-1] == "hdf") and (
                    MOD11_file.split(".")[2] in [tile_Key]
                ):
                    hdf_list.append(join(data_path, MOD11Folder, MOD11_file))
        else:
            continue

    return hdf_list


"""
function to convert date(YYYYMMDD) to DOY
"""


def date2DOY(file_date):
    file_date = datetime.strptime(file_date, "%Y%m%d")
    days = str(file_date.year).zfill(4) + str(file_date.timetuple().tm_yday).zfill(3)
    return days


"""
function to extract MODIS LST data given latitude, longitude and a period of time
input of the function:
points: coordinate (lat, lon); float, a list of tuple
startDate: start time, should be in the format of 'YYYYMMDD' [string]
endDate: end time, should be in the format of "YYYYMMDD" [string]
output of the function:
a list of data with format[{coordinate: value, 'MOD11A1': { day/night: [(time1,value1),
(time2,value2),...]}, 'MYD11A1': { day/night: [(time1,value1),(time2,value2),...]}},...]
---------------Valid LST Data-----------------
7500 <= LST <= 65535 (scale is 0.02)
QC == 0 or QC & 0x000F == 1
--------------Invalid LST Data----------------
LST = -1
Others(eg. QC == 2/3)
"""


def extract_MODIS_LST(tile_dict, startDate, endDate):
    # Results Initialization
    MODIS_LST_res = {
        ("{:.6f}".format(pt[0]), "{:.6f}".format(pt[1])): {"MOD11A1": [], "MYD11A1": []}
        for pt in points
    }

    # MODIS LST folders
    MOD_data_path = join(os.path.expanduser("~"), settings.MOD11A1_PATH)
    MYD_data_path = join(os.path.expanduser("~"), settings.MYD11A1_PATH)

    # Get Data
    for tile_Key in tile_dict.keys():
        hdf_files = []
        hdf_files = hdf_files + create_tar_hdf(
            MOD_data_path, tile_Key, startDate, endDate
        )
        hdf_files = hdf_files + create_tar_hdf(
            MYD_data_path, tile_Key, startDate, endDate
        )

        for file_path in hdf_files:
            try:
                hdf_ds = gdal.Open(file_path, gdal.GA_ReadOnly)
            except Exception:
                print("Unable to open file: \n" + file_path + "\nSkipping\n")
                continue

            # Compare file info and folder info
            path_info = file_path.split("/")
            file_date = path_info[-2].replace(".", "")
            file_type = path_info[-3].split(".")[0]
            file_name = os.path.basename(file_path)
            file_info = file_name.split(".")
            if file_type != file_info[0] or date2DOY(file_date) != file_info[1][1:]:
                print("File info unmatch: \n" + file_path + "\nSkipping\n")
                continue

            # read dataset
            # Get Day LST and QC
            # Get Day QC Dataset
            day_lst_ds = gdal.Open(hdf_ds.GetSubDatasets()[0][0])
            day_lst_rst = day_lst_ds.GetRasterBand(1).ReadAsArray()
            day_qc_ds = gdal.Open(hdf_ds.GetSubDatasets()[1][0])
            day_qc_rst = day_qc_ds.GetRasterBand(1).ReadAsArray()

            # Get Night LST and QC
            # Get Night QC Value
            night_lst_ds = gdal.Open(hdf_ds.GetSubDatasets()[4][0])
            night_lst_rst = night_lst_ds.GetRasterBand(1).ReadAsArray()
            night_qc_ds = gdal.Open(hdf_ds.GetSubDatasets()[5][0])
            night_qc_rst = night_qc_ds.GetRasterBand(1).ReadAsArray()

            for location in tile_dict[tile_Key]:
                latitude = location[0]
                longitude = location[1]

                # Get img coords(from proj coord) and check if data is valid
                is_valid, px, py = geo_functions.point_boundary_is_valid(
                    day_qc_ds, latitude, longitude
                )
                if not is_valid:
                    continue

                # LST value result
                valid_lst_value = []

                # Get Day QC Value
                day_lst_value = -1
                day_qc_value = geo_functions.get_band_value(day_qc_rst, px, py, 3, True)
                if day_qc_value == 0 or day_qc_value & 0x000F == 1:
                    # Get Day LST Value
                    day_lst_value = geo_functions.get_band_value(day_lst_rst, px, py, 3)
                # Add to valid list
                if 7500 <= day_lst_value <= 65535:
                    valid_lst_value.append(day_lst_value)

                # Get Night QC Value
                night_lst_value = -1
                night_qc_value = geo_functions.get_band_value(
                    night_qc_rst, px, py, 3, True
                )
                if night_qc_value == 0 or night_qc_value & 0x000F == 1:
                    # Get Night LST Value
                    night_lst_value = geo_functions.get_band_value(
                        night_lst_rst, px, py, 3
                    )
                # Add to valid list
                if 7500 <= night_lst_value <= 65535:
                    valid_lst_value.append(night_lst_value)

                # Get Final LST Value
                lst_value = -1
                if len(valid_lst_value) != 0:
                    lst_value = np.mean(np.array(valid_lst_value))

                # Add to results
                res_key = geo_functions.get_latlon_key(latitude, longitude)

                MODIS_LST_record = MODIS_LST_res[res_key]
                MODIS_LST_record[file_type].append(
                    (file_date, float("%.2f" % lst_value))
                )

    # sort the results by Date
    for _, MODIS_LST_record in MODIS_LST_res.items():
        if "MOD11A1" in MODIS_LST_record.keys():
            MODIS_LST_record["MOD11A1"].sort()
        if "MYD11A1" in MODIS_LST_record.keys():
            MODIS_LST_record["MYD11A1"].sort()

    return MODIS_LST_res


# test main function
if __name__ == "__main__":
    # points = [(47.729849, -119.540372), (43.79562, -111.466202),
    # (42.034932, -107.947733), (30.732852, -104.768943), (42.969349, -109.885573)]
    # points = list(np.load('lonlatArray.npz')['arr_0'])
    points = (((np.load("2016TT_points.npz")["arr_0"]).tolist())["Corn"])[0:]
    # points = [(36.822911628984635, -90.42548077091972)]
    startDate = "20160401"
    endDate = "20161001"
    MODIS_LST_res = extract_MODIS_LST(points, startDate, endDate)
    printer.pprint(MODIS_LST_res)
    print("over")
