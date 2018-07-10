# coding: utf8
import os, sys, pprint
from osgeo import gdal,osr
from landsat_index import ConvertToWRS
from os.path import join
import glob, time, struct
from datetime import date
from gdalconst import *
import settings
import numpy as np

sys.path.append(join(os.path.dirname(os.path.realpath(__file__)),".."))
printer = pprint.PrettyPrinter(indent=3)

band_key_list = ['R_band', 'G_band', 'B_band', 'NIR_band', 'SWIR1', 'SWIR2']

'''
calculate the path and row of each point
'''
def target_path_row(sample_list):
    pr = list()
    result_sample = {}
    to_WRS = ConvertToWRS()
    for (lat,lon) in (sample_list):
        result = to_WRS.get_wrs(lat,lon)
        for i in range(0,len(result)):
            path = '0' + str(result[i]['path'])
            row = '0' + str(result[i]['row'])
            pr.append((path,row))
            key = path + '-' + row
            result_sample.setdefault(key,[]).append((lat,lon))
    return pr, result_sample

'''
get the suitable file folds acoording to the path and row
'''
def get_file_list(path_row, data_dir, start_time, end_time):
    tar_list = {}
    for (path,row) in path_row:
        file_folds = glob.glob(join(data_dir,"*","01",path,row,"*"))

        file_list = []
        for file_folder in file_folds:
            file_date = file_folder.split("_")[-4]
            if start_time <= file_date <= end_time:
                file_list.append(file_folder)

        key = path + "-" + row
        tar_list[key] = (file_list)
    return tar_list

'''
the data format supported by gdal
'''
def pt2fmt(pt):
    fmttypes = {GDT_Byte: 'B',GDT_Int16: 'h',GDT_UInt16: 'H',GDT_Int32: 'i',\
    GDT_UInt32: 'I',GDT_Float32: 'f',GDT_Float64: 'f'}
    return fmttypes.get(pt, 'x')

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

def lonlat2geo(dataset, lon, lat):
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


"""
function to extract BLUE/GREEN/RED/NIR/SWIR1/SWIR2 data given sample points and a period of time
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
        value = round(result*0.0001,4)
    else:
        value = round(result,4)
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

def extract_Landsat_SR(points, start_time, end_time):
    time0 = time.time()
    path_row , sample_points = target_path_row(points)
    if len(path_row) == 0:
        raise Exception("There is no consitent path and row")
    else:
        pass

    path_row = list(set(path_row))
    print(path_row)
    time1 = time.time()
    print('Get Path_Row time: ' + str(time1 - time0))

    #Results Initialization
    landsat_ref_res = []
    for pt in points:
        landsat_ref_record = {}
        #Add coordinate key
        landsat_ref_record['coordinate'] = (pt[0], pt[1])
        landsat_ref_record['B_band'] = []
        landsat_ref_record['G_band'] = []
        landsat_ref_record['R_band'] = []
        landsat_ref_record['NIR_band'] = []
        landsat_ref_record['SWIR1'] = []
        landsat_ref_record['SWIR2'] = []
        landsat_ref_res.append(landsat_ref_record)

    tar_list = {}
    for landsat_dir in settings.LANDSAT_SR_PATH:
        data_dir = join(os.path.expanduser('~'), landsat_dir)
        tar_mid = get_file_list(path_row, data_dir, start_time, end_time)
        for k,v in tar_mid.items():
            if k in tar_list.keys():
                tar_list[k] += v
            else:
                tar_list[k] = v

    if len(tar_list) == 0:
        raise Exception("There is no suitable files")
    else:
        pass

    time2 = time.time()
    print('Initialization time: ' + str(time2 - time1))

    #Try to extract Data in a faster way
    for pr_key in tar_list.keys():
        img_folders = tar_list[pr_key]
        for folder_path in img_folders:
            file_date = folder_path.split("_")[-4]
            #open qc img
            try:
                qc_path = join(folder_path,folder_path.split("/")[-1] + '_pixel_qa.img')
                msband = gdal.Open(qc_path)
            except:
                print('Unable to open QC file in folder:' + folder_path + '\n')
                continue

            maskdata = msband.GetRasterBand(1)
            geomatrix = msband.GetGeoTransform()
            xsize = msband.RasterXSize
            ysize = msband.RasterYSize

            time3 = time.time()
            for location in sample_points[pr_key]:
                lat = location[0]
                lon = location[1]
                x, y = lonlat2geo(msband, lon, lat)
                px = int((x - geomatrix[0])/geomatrix[1])
                py = int((y - geomatrix[3])/geomatrix[5])
                if px > xsize or px < 0 or py > ysize or py < 0:
                    printer.pprint("the point is out of the image range(" + folder_path + "):" + str(location))
                    continue
                else:
                    pass

                B_ref = -1.0000
                G_ref = -1.0000
                R_ref = -1.0000
                N_ref = -1.0000
                SW1_ref = -1.0000
                SW2_ref = -1.0000

                cloud_mask = maskdata.ReadRaster(px, py, 1, 1, buf_type = maskdata.DataType)
                fmt2 = pt2fmt(maskdata.DataType)
                intval2 = struct.unpack(fmt2, cloud_mask)

                time31 = time.time()
                # print('Find Cloudmask time:' + str(time31 - time3))

                #For clear sky, open other img files
                if intval2[0] in [66, 130] or intval2[0] in [322, 386, 834, 898, 1346]:
                    #Open band files
                    try:
                        satellite = folder_path.split("/")[-5]
                        if satellite in ["LT05","LE07"]:
                            B_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band1.img')
                            G_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band2.img')
                            R_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band3.img')
                            N_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band4.img')
                            SW1_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band5.img')
                            SW2_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band7.img')
                        else:
                            B_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band2.img')
                            G_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band3.img')
                            R_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band4.img')
                            N_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band5.img')
                            SW1_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band6.img')
                            SW2_path = join(folder_path,folder_path.split("/")[-1] + '_sr_band7.img')

                        #Get B Band Ref
                        bandfile = gdal.Open(B_path)
                        banddata = bandfile.GetRasterBand(1)
                        result = banddata.ReadRaster(px, py, 1, 1, buf_type=banddata.DataType)
                        fmt = pt2fmt(banddata.DataType)
                        intval = struct.unpack(fmt, result)
                        B_ref = round(intval[0]*0.0001,4)

                        #Get G Band Ref
                        bandfile = gdal.Open(G_path)
                        banddata = bandfile.GetRasterBand(1)
                        result = banddata.ReadRaster(px, py, 1, 1, buf_type=banddata.DataType)
                        fmt = pt2fmt(banddata.DataType)
                        intval = struct.unpack(fmt, result)
                        G_ref = round(intval[0]*0.0001,4)

                        #Get R Band Ref
                        bandfile = gdal.Open(R_path)
                        banddata = bandfile.GetRasterBand(1)
                        result = banddata.ReadRaster(px, py, 1, 1, buf_type=banddata.DataType)
                        fmt = pt2fmt(banddata.DataType)
                        intval = struct.unpack(fmt, result)
                        R_ref = round(intval[0]*0.0001,4)

                        #Get N Band Ref
                        bandfile = gdal.Open(N_path)
                        banddata = bandfile.GetRasterBand(1)
                        result = banddata.ReadRaster(px, py, 1, 1, buf_type=banddata.DataType)
                        fmt = pt2fmt(banddata.DataType)
                        intval = struct.unpack(fmt, result)
                        N_ref = round(intval[0]*0.0001,4)

                        #Get SW1 Band Ref
                        bandfile = gdal.Open(SW1_path)
                        banddata = bandfile.GetRasterBand(1)
                        result = banddata.ReadRaster(px, py, 1, 1, buf_type=banddata.DataType)
                        fmt = pt2fmt(banddata.DataType)
                        intval = struct.unpack(fmt, result)
                        SW1_ref = round(intval[0]*0.0001,4)

                        #Get SW2 Band Ref
                        bandfile = gdal.Open(SW2_path)
                        banddata = bandfile.GetRasterBand(1)
                        result = banddata.ReadRaster(px, py, 1, 1, buf_type=banddata.DataType)
                        fmt = pt2fmt(banddata.DataType)
                        intval = struct.unpack(fmt, result)
                        SW2_ref = round(intval[0]*0.0001,4)

                        # print('Extractor success!' + folder_path + '\n')

                    except:
                        print("Cannot open img files in folder:" + folder_path + "\n")
                        continue


                time32 = time.time()
                print('Extract data time:' + str(time32 - time31))
                #Add to results
                point_exist = False
                for landsat_ref_record in landsat_ref_res:
                    if landsat_ref_record['coordinate'] == location:
                        point_exist = True
                        try:
                            landsat_ref_record['B_band'].append((file_date, float(B_ref)))
                            landsat_ref_record['G_band'].append((file_date, float(G_ref)))
                            landsat_ref_record['R_band'].append((file_date, float(R_ref)))
                            landsat_ref_record['NIR_band'].append((file_date, float(N_ref)))
                            landsat_ref_record['SWIR1'].append((file_date, float(SW1_ref)))
                            landsat_ref_record['SWIR2'].append((file_date, float(SW2_ref)))
                            break
                        except KeyError:
                            pass

                if point_exist == False:
                    print("The point has not been initialized: " + str(location) + "\n")
                    continue

                time33 = time.time()
                # print('Add result time:' + str(time33 - time32))

            time4 = time.time()
            print('one img time:' + str(time4 - time3))

    time5 = time.time()
    #Sort the resuls by date
    for landsat_ref_record in landsat_ref_res:
        landsat_ref_record['B_band'].sort()
        landsat_ref_record['G_band'].sort()
        landsat_ref_record['R_band'].sort()
        landsat_ref_record['NIR_band'].sort()
        landsat_ref_record['SWIR1'].sort()
        landsat_ref_record['SWIR2'].sort()

    time6 = time.time()
    print('Sort result time:' + str(time6 - time5))

    return landsat_ref_res

if __name__ == '__main__':

    startt=time.time()
    points = [(29.822911628984635, -100.42548077091972)]
    start_time = "20130301"
    end_time = "20131030"
    result = extract_Landsat_SR(points, start_time, end_time)
    printer.pprint(result)
    print(time.time() - startt)
