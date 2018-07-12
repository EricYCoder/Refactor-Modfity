import os, sys, pprint
import json, glob, time, struct
import numpy as np
from os.path import join
from datetime import date
from osgeo import gdal, osr
from sentinel_index import ConvertToMRGS
from gdalconst import *


printer = pprint.PrettyPrinter(indent=3)
home_dir = os.path.expanduser('~')

# load process sentinel list, it must be update 
file_json = join(os.path.expanduser('~'), 'waterfall/total_list.json')
with open(file_json, 'r') as fp:
    all_file_folds = json.load(fp)

'''
calculate the tiler of each point
'''


def get_sentinel_tile(points):
    tiles = list()
    tile_dict = {}
    to_mgrs = ConvertToMRGS()
    for (lat, lon) in points:
        result = to_mgrs.get_mgrs(lat, lon)
        for key in result:
            tiles.append(key)
            tile_dict.setdefault(key, []).append((lat, lon))
    return tiles, tile_dict


'''
get the suitable file folds acoording to the MGRS
'''


def get_file_list(tiles, start_time, end_time):
    tar_list = {}
    count = 0
    if len(tiles) > 0:
        for tile in tiles:
            tmp_tile = join(tile[0:2], tile[2:3], tile[3:5])
            file_folds = [(home_dir + tmp) for tmp in all_file_folds if tmp_tile in tmp]  # find the file in list
            file_list = []
            for file_fold in file_folds:
                file_date = file_fold.split('/')[9] + '%02d' % int(file_fold.split('/')[10]) + '%02d' % int(file_fold.split('/')[11])
                if start_time <= file_date <= end_time:
                    file_fold = glob.glob(join(file_fold, 'GRANULE', '*', 'IMG_DATA'))[0]
                    file_list.append(file_fold)
            tar_list[tile] = (file_list)
            count += 1
            print(count/len(tiles))
    else:
        print('Index is wrong and please check!')

    return tar_list

def getSRSPair(dataset):
    '''
    获得给定数据的投影参考系和地理参考系
    :param dataset: GDAL地理数据
    :return: 投影参考系和地理参考系
    '''
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs

def latlon2geo(dataset, lat, lon):
    '''
    将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param lon: 地理坐标lon经度
    :param lat: 地理坐标lat纬度
    :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lon, lat)
    return coords[:2]

def geo2imagexy(dataset, x, y):
    '''
    根据GDAL的六参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
    :param dataset: GDAL地理数据
    :param x: 投影或地理坐标x
    :param y: 投影或地理坐标y
    :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col)
    '''
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b) # 使用numpy的linalg.solve进行二元一次方程的求解

"""
function to extract data given sample points and a period of time
input of the function:
points: coordinate (lat, lon); float, a list of tuple
startDate: start time, should be in the format of 'YYYYMMDD' [string]
endDate: end time, should be in the format of "YYYYMMDD" [string]
output of the function:
a list of data with format[{'coordinate': (lat,lon),'BLUE': [(time1,value1),(time2,value2),....], 'GREEN': [(time1,value1),(time2,value2),....],...}...]
--------------valid Data----------------
0-1
--------------Invalid Data----------------
-1
"""

def get_band_value(bandraster, px, py, cloud=False):
    result = bandraster[py][px]
    if cloud == False:
        value = round(result*0.0001, 4)
    else:
        value = round(result, 4)
    return value

def point_boundary_is_valid(bandfile, lat, lon):
    xsize = bandfile.RasterXSize
    ysize = bandfile.RasterYSize
    x, y = latlon2geo(bandfile, lat, lon) # (lat, lon) to (x, y)
    px, py = geo2imagexy(bandfile, x, y)  # (x, y) to (row, col)
    px = int(px)
    py = int(py)
    if px > xsize or px < 0 or py > ysize or py < 0:
        printer.pprint("the point is out of the image range:" + str(location))
        return False, px, py
    else:
        return True, px, py

def extract_sentinel_SR(points, start_time, end_time):
    global_count = 0

    tiles, sample_points = get_sentinel_tile(points)
    if len(tiles) == 0:
        raise Exception("There is no consitent path and row")
    else:
        pass

    # Results Initialization
    sentinel_ref_res = {('{:.6f}'.format(pt[0]), '{:.6f}'.format(pt[1]): {"B_band":[],"G_band":[],
                        "R_band":[],"NIR_band":[],"SWIR1":[],"SWIR2":[]} for pt in points}

    # get tiles list
    print('Getting file list')
    tar_list = get_file_list(tiles, start_time, end_time)

    if len(tar_list) == 0:
        raise Exception("There is no suitable files")
    else:
        pass

    tile_count = 0
    for pr_key in tar_list.keys():
        
        # display the message
        time_stamp = datetime.datetime.now()
        print('System time:', time_stamp.strftime("%Y-%m-%d %H:%M:%S %p")
        print('Processing the tile', tile_count, 'out of ', len(tar_list.keys()), ', file number:' , len(tar_list[pr_key]))
   
        tile_count += 1

        img_folders = tar_list[pr_key]
        
        for folder_path in img_folders:
            file_date = folder_path.split('/')[9]+ '%02d' % int(folder_path.split('/')[10]) + '%02d' % int(folder_path.split('/')[11])
     
            try:
                # creat cloud path
                qc_path = folder_path.split('/S2')[0].replace('SAFE_sentinel/','SAFE_sentinel/cloudmask/sentinel/')
                if os.path.exists(qc_path):
                    if len(list(os.listdir(qc_path))) != 4:
                        continue
                else:
                    continue

                try:
                    cloudmask = gdal.Open(join(qc_path, 'cloud.img'))
                    cloudmask_raster = cloudmask.GetRasterBand(1).ReadAsArray()
                except Exception as e:
                    print(e)
                    print('Unable to open QC file in folder:' + qc_path + '\n')
                    continue

            except Exception as e:
                print(e)
                print('Unable again to open QC file in folder:' + folder_path + '\n')
                continue

            # get all data path for 20m
            R20_list = os.listdir(join(folder_path, 'R20m'))
            path_list = R20_list
            bandfiles = {}
            bandraster = {}

            if len(path_list) > 0: # [0]set the data path
                for tmp_path in path_list:
                    try:
                        if '_B02_20m' in tmp_path:
                            b2_path = join(folder_path, 'R20m', tmp_path)
                            bandfiles['B_band'] = gdal.Open(b2_path)
                            bandraster['B_band'] = bandfiles['B_band'].GetRasterBand(1).ReadAsArray()
                        elif '_B03_20m' in tmp_path:
                            b3_path = join(folder_path, 'R20m', tmp_path)
                            bandfiles['G_band'] = gdal.Open(b3_path)
                            bandraster['G_band'] = bandfiles['G_band'].GetRasterBand(1).ReadAsArray()
                        elif '_B04_20m' in tmp_path:
                            b4_path = join(folder_path, 'R20m', tmp_path)
                            bandfiles['R_band'] = gdal.Open(b4_path)
                            bandraster['R_band'] = bandfiles['R_band'].GetRasterBand(1).ReadAsArray()
                        elif '_B8A_20m' in tmp_path:
                            b8a_path = join(folder_path, 'R20m', tmp_path)
                            bandfiles['NIR_band'] = gdal.Open(b8a_path)
                            bandraster['NIR_band'] = bandfiles['NIR_band'].GetRasterBand(1).ReadAsArray()
                        elif '_B11_20m' in tmp_path:
                            b11_path = join(folder_path, 'R20m', tmp_path)
                            bandfiles['SWIR1'] = gdal.Open(b11_path)
                            bandraster['SWIR1'] = bandfiles['SWIR1'].GetRasterBand(1).ReadAsArray()
                        elif '_B12_20m' in tmp_path:
                            b12_path = join(folder_path, 'R20m', tmp_path)
                            bandfiles['SWIR2'] = gdal.Open(b12_path)
                            bandraster['SWIR2'] = bandfiles['SWIR2'].GetRasterBand(1).ReadAsArray()
                        else:
                            continue
                    except AttributeError as e:
                        print(e)
                        continue
            else:
                print("Process fail!")
                continue

            for location in sample_points[pr_key]:
                # print(location)
                lat = location[0]
                lon = location[1]
                band_data = {
                    'B_band': -1,
                    'R_band': -1,
                    'G_band': -1,
                    'NIR_band': -1,
                    'SWIR1': -1,
                    'SWIR2': -1
                }

                is_valid, px, py = point_boundary_is_valid(cloudmask, lat, lon)
                if not is_valid:
                    continue
                cloudmask_data = get_band_value(cloudmask_raster, px, py, cloud=True)
                if cloudmask_data !=1:
                    continue

                # get 20m data
                try:
                    bandfile = bandfiles['NIR_band']
                    is_valid, px, py = point_boundary_is_valid(bandfile, lat, lon)
                    if not is_valid:
                        continue
                    band_data['NIR_band'] = get_band_value(bandraster['NIR_band'], px, py)
                    band_data['SWIR1'] = get_band_value(bandraster['SWIR1'], px, py)
                    band_data['SWIR2'] = get_band_value(bandraster['SWIR2'], px, py)
                    band_data['B_band'] = get_band_value(bandraster['B_band'], px, py)
                    band_data['G_band'] = get_band_value(bandraster['G_band'], px, py)
                    band_data['R_band'] = get_band_value(bandraster['R_band'], px, py)
                except:
                    print("Cannot open 20m img files in folder:" + folder_path + "\n")


                # #Add to results
                # point_exist = False
                # for sentinel_ref_record in sentinel_ref_res[crop_type]:
                #     if sentinel_ref_record['coordinate'] == location:
                #         point_exist = True
                #         for band_type in band_key_list:
                #             try:
                #                 sentinel_ref_record[band_type].append((file_date, float(band_data[band_type])))
                #                 global_count += 1
                #                 if global_count % 3000:
                #                     global_toc = time.time()
                #                     # print('Speed of data per day per band is', (global_toc-global_tic)/global_count)
                #             except KeyError:
                #                 print('Key Error!!!')
                #                 pass
                #         break
                # if point_exist == False:
                #     # print("The point has not been initialized: " + str(location))
                #     continue

    #Sort the resuls by date
    for sentinel_ref_record in sentinel_ref_res[crop_type]:
        for band_type in band_key_list:
            sentinel_ref_record[band_type].sort()

    return sentinel_ref_res

if __name__ == '__main__':
    startt=time.time()
    point_json_corn = join(os.path.expanduser('~'), 
                    'data_pool/waterfall_data/2017_Corn_points.json')
    point_json_soybeans = join(os.path.expanduser('~'),
                    'data_pool/waterfall_data/2017_Soybeans_points.json')
    point_json_other = join(os.path.expanduser('~'),
                    'data_pool/waterfall_data/2017_Other_points.json')
    with open(point_json_corn , 'r') as fp:
        corn_points = json.load(fp)
    with open(point_json_soybeans , 'r') as fp:
        soybeans_points = json.load(fp)
    with open(point_json_other, 'r') as fp:
        other_points = json.load(fp)

    points = {
        'Corn': corn_points[:2],
        'Soybeans': soybeans_points[:2],
        'Other': other_points[:2]
    }
    start_time = "20170501"
    end_time = "20170527"

    result = extract_sentinel_SR(points, start_time, end_time)

    for crop_type in ['Corn', 'Soybeans', 'Other']:
        np.savez('2017_Corn_Sentinel.npz', result['Corn'])
        np.savez('2017_Soybeans_Sentinel.npz', result['Soybeans'])
        np.savez('2017_Other_Sentinel.npz', result['Other'])
    print(time.time() - startt)
