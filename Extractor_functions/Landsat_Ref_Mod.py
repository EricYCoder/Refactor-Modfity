# coding: utf8
import os
import sys
import pprint
from osgeo import gdal
from os.path import join
import glob
import settings
import geo_functions
import numpy as np

sys.path.append(join(os.path.dirname(os.path.realpath(__file__)), ".."))
printer = pprint.PrettyPrinter(indent=3)

bandNum = len(settings.band_key_list)
LT57_band_dict = {
    settings.band_key_list[i]: settings.LT57_band_index[i] for i in range(bandNum)
}
LT8_band_dict = {
    settings.band_key_list[i]: settings.LT8_band_index[i] for i in range(bandNum)
}

"""
get the suitable file folds acoording to the path and row
"""


def get_file_list(path_rows, data_dirs, start_time, end_time):
    tar_list = {}
    for (path, row) in path_rows:
        for data_dir in data_dirs:
            file_folds = glob.glob(os.path.join(data_dir, "*", "01", path, row, "*"))

            file_list = []
            for file_folder in file_folds:
                file_date = file_folder.split("_")[-4]
                if start_time <= file_date <= end_time:
                    file_list.append(file_folder)

            key = path + "-" + row
            if key in tar_list.keys():
                tar_list[key] += file_list
            else:
                tar_list[key] = file_list
    return tar_list


"""
function to extract BLUE/GREEN/RED/NIR/SWIR1/SWIR2 data given sample points
and a period of time
input of the function:
points: coordinate (lat, lon); float, a list of tuple
startDate: start time, should be in the format of 'YYYYMMDD' [string]
endDate: end time, should be in the format of "YYYYMMDD" [string]
output of the function:
a list of data with format[{'coordinate': (lat,lon),'BLUE': [(time1,value1),
(time2,value2),....], 'GREEN': [(time1,value1),(time2,value2),....],...}...]
--------------valid Data----------------
0-1
--------------Invalid Data----------------
-1
"""


def extract_Landsat_SR(tile_dict, start_time, end_time):
    # Get path_rows in tile_dict
    path_rows = tile_dict.keys()
    # Results
    landsat_ref_res = {}

    tar_list = {}
    data_dirs = []
    for landsat_dir in settings.LANDSAT_SR_PATH:
        data_dirs.append(join(os.path.expanduser("~"), landsat_dir))
    tar_list = get_file_list(path_rows, data_dirs, start_time, end_time)

    if len(tar_list) == 0:
        raise Exception("There is no suitable files")
        return
    else:
        pass

    tile_count = 0
    for pr_key in tar_list.keys():
        print("Processing the tile", tile_count, "out of ", len(tar_list.keys()))
        tile_count += 1

        img_folders = tar_list[pr_key]

        for folder_path in img_folders:
            file_date = folder_path.split("_")[-4]

            # open qc img
            try:
                qc_path = join(
                    folder_path, folder_path.split("/")[-1] + "_pixel_qa.img"
                )
                cloudmask = gdal.Open(qc_path)
                cloudmask_raster = cloudmask.GetRasterBand(1).ReadAsArray()
            except Exception as e:
                print(e)
                print("Unable to open QC file\n")
                continue

            # open band img
            satellite = folder_path.split("/")[-5]
            File_Path = {}
            bandfiles = {}
            bandraster = {}

            if satellite in ["LT05", "LE07"]:
                for band_type in settings.band_key_list:
                    File_Path[band_type] = join(
                        folder_path,
                        folder_path.split("/")[-1] + LT57_band_dict[band_type],
                    )
            elif satellite in ["LC08"]:
                for band_type in settings.band_key_list:
                    File_Path[band_type] = join(
                        folder_path,
                        folder_path.split("/")[-1] + LT8_band_dict[band_type],
                    )
            else:
                print("satellite type error!")
                continue

            for band_type in settings.band_key_list:
                try:
                    bandfiles[band_type] = gdal.Open(File_Path[band_type])
                    bandraster[band_type] = (
                        bandfiles[band_type].GetRasterBand(1).ReadAsArray()
                    )

                except Exception as e:
                    print(e)
                    print("Unable to open " + band_type + " file\n")
                    continue

            for location in tile_dict[pr_key]:
                lat = location[0]
                lon = location[1]

                # Check if data is valid
                is_valid, px, py = geo_functions.point_boundary_is_valid(
                    cloudmask, lat, lon
                )
                if not is_valid:
                    continue
                cloudmask_data = geo_functions.get_band_value(
                    cloudmask_raster, px, py, 1, cloud=True
                )
                if cloudmask_data not in [66, 130, 322, 386, 834, 898, 1346]:
                    continue

                band_data = {band_type: -1 for band_type in settings.band_key_list}
                # Get Band Ref
                for band_type in settings.band_key_list:
                    try:
                        ref_value = geo_functions.get_band_value(
                            bandraster[band_type], px, py, 1
                        )
                        if 0.0 < ref_value < 1.0:
                            band_data[band_type] = ref_value
                        else:
                            band_data[band_type] = -1.0
                    except Exception:
                        print("Unable to get " + band_type + " data\n")
                        continue

                # Add to results
                res_key = geo_functions.get_latlon_key(lat, lon)
                if res_key in landsat_ref_res.keys:
                    landsat_ref_record = landsat_ref_res[res_key]
                    for band_type in settings.band_key_list:
                        try:
                            landsat_ref_record[band_type].append(
                                (file_date, float(band_data[band_type]))
                            )
                        except KeyError:
                            print("Key Error!!!")
                            pass
                else:
                    landsat_ref_record = {}
                    for band_type in settings.band_key_list:
                        landsat_ref_record[band_type] = (
                            file_date,
                            float(band_data[band_type]),
                        )
                    landsat_ref_res[res_key] = landsat_ref_record

    # Sort the resuls by date
    for _, landsat_ref_record in landsat_ref_res.items():
        for band_type in settings.band_key_list:
            landsat_ref_record[band_type].sort()

    return landsat_ref_res


if __name__ == "__main__":
    print("begin!")
    # points = [(40.570, -89.354), (42.353, -89.352), (40.619, -116.959),
    #           (48.687, -96.205), (41.277, -95.359), (42.375, -83.369),
    #          (36.968, -77.950), (34.65, -75.35), (41.46, -84.53), (25.25, -77.95),
    #           (35.22, -71.12), (40.565, -77.359),
    #          (42.347, -95.351), (43.290, -98.182), (37.901, -99.482),
    #           (40.569, -95.354), (42.361, -96.506), (42.361, -96.606)]
    # points = [(36.822911628984635, -90.42548077091972)]
    # points = np.load('lonlatArray.npz')['arr_0']
    points = (((np.load("2016TT_points.npz")["arr_0"]).tolist())["Corn"])[0:]

    start_time = "20160401"
    end_time = "20161001"
    result = extract_Landsat_SR(points, start_time, end_time)
    printer.pprint(result)
    print("finish!")
