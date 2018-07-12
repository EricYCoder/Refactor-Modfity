from osgeo import osr
import numpy as np

"""
获得给定数据的投影参考系和地理参考系
:param dataset: GDAL地理数据
:return: 投影参考系和地理参考系
"""


def getSRSPair(dataset):
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs


"""
将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
:param dataset: GDAL地理数据
:param lon: 地理坐标lon经度
:param lat: 地理坐标lat纬度
:return: 经纬度坐标(lon, lat)对应的投影坐标
"""


def lonlat2geo(dataset, lon, lat):
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lon, lat)
    return coords[:2]


"""
根据GDAL的六参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
:param dataset: GDAL地理数据
:param x: 投影或地理坐标x
:param y: 投影或地理坐标y
:return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col)
"""


def geo2imagexy(dataset, x, y):
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b)  # 使用numpy的linalg.solve进行二元一次方程的求解


"""
根据图像坐标获取波段内像素值
:param bandraster: 栅格图像波段
:param px: 图像坐标x
:param py: 图像坐标y
:param dataType: 获取数据类型（1:Landsat 2:Sentinel 3:MODIS）
:param cloud: 是否是云掩膜
"""


def get_band_value(bandraster, px, py, dataType, cloud=False):
    result = bandraster[py][px]
    value = -1
    if cloud:
        value = round(result, 4)
    elif dataType is 1 or 2:
        value = round(result * 0.0001, 4)
    elif dataType is 3:
        value = round(result * 0.02, 4)
    else:
        value = -1
    return value


"""
根据经纬度获取图像坐标
:param bandfile: 栅格图像文件句柄
:param lat: 纬度
:param lon: 经度
"""


def point_boundary_is_valid(bandfile, lat, lon):
    xsize = bandfile.RasterXSize
    ysize = bandfile.RasterYSize
    x, y = lonlat2geo(bandfile, lon, lat)  # (lat, lon) to (x, y)
    px, py = geo2imagexy(bandfile, x, y)  # (x, y) to (row, col)
    px = int(px)
    py = int(py)
    if px > xsize or px < 0 or py > ysize or py < 0:
        print("the point is out of the image range:" + str(lat) + ", " + str(lon))
        return False, px, py
    else:
        return True, px, py


"""
根据经纬度获取结果字典key
:param lat: 纬度
:param lon: 经度
"""


def get_latlon_key(lat, lon):
    res_key = ("{:.6f}".format(lat), "{:.6f}".format(lon))
    return res_key
