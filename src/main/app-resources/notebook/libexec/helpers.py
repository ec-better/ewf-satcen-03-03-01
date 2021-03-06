from py_snap_helpers import *
import otbApplication
from ogr import osr
import ogr
import gdal
import geopandas as gp
import os
import pandas as pd
import numpy as np


from shapely.geometry import box
from shapely.wkt import loads
import sys

gdal.UseExceptions()

#sys.path.append('/opt/OTB/lib/python')
#sys.path.append('/opt/OTB/lib/libfftw3.so.3')
#sys.path.append('/opt/anaconda/bin')
#os.environ['OTB_APPLICATION_PATH'] = '/opt/OTB/lib/otb/applications'
#os.environ['LD_LIBRARY_PATH'] = '/opt/OTB/lib'
#os.environ['ITK_AUTOLOAD_PATH'] = '/opt/OTB/lib/otb/applications'
os.environ['PREFIX'] = '/opt/anaconda/envs/env_ewf_satcen_03_03_01/'
sys.path.append(os.path.join(os.environ['PREFIX'], 'bin'))

os.environ['_JAVA_OPTIONS'] = '-Xms24g -Xmx24g'

from gdal_calc import Calc as gdalCalc

import cioppy
ciop = cioppy.Cioppy()


def check_existence(product):

    if os.path.isfile('{}'.format(product)):
        print('output file: {} generated successfully.'.format(product))
    else:
        raise RuntimeError('Missing output file: {}.dim'.format(product))



def get_metadata(input_references, data_path):

    if isinstance(input_references, str):

        search_params = dict()

        search_params['do'] = 'terradue'

        products = gp.GeoDataFrame(ciop.search(end_point=input_references, 
                                               params=search_params,
                                               output_fields='identifier,self,wkt,startdate,enddate,enclosure,orbitDirection,track,orbitNumber',
                                               timeout=360000,
                                               model='EOP'))

    else:    

        temp_results = []

        for index, self in enumerate(input_references):

            search_params = dict()

            search_params['do'] = 'terradue'

            temp_results.append(ciop.search(end_point=self, 
                                params=search_params,
                                output_fields='identifier,self,wkt,startdate,enddate,enclosure,orbitDirection,track,orbitNumber', 
                                model='EOP')[0])

        products = gp.GeoDataFrame(temp_results)
        
        products = products.merge(products.apply(lambda row: analyse(row, data_path), axis=1),
                                    left_index=True,
                                  right_index=True)
        
        
    return products

def analyse(row, data_path):
    
    series = dict()

    series['local_path'] = os.path.join(data_path, row['identifier'], row['identifier'] + '.SAFE', 'manifest.safe')
   
    return pd.Series(series)

def bbox_to_wkt(bbox):
    
    return box(*[float(c) for c in bbox.split(',')]).wkt


def get_inteserction_aoi_prod(aoi,prod_wkt):
    
    return loads(prod_wkt).intersection(loads(aoi)).wkt


def zone(coordinates):
    
    if 56 <= coordinates[1] < 64 and 3 <= coordinates[0] < 12:
        return 32
    if 72 <= coordinates[1] < 84 and 0 <= coordinates[0] < 42:
        if coordinates[0] < 9:
            return 31
        elif coordinates[0] < 21:
            return 33
        elif coordinates[0] < 33:
            return 35
        return 37
    return int((coordinates[0] + 180) / 6) + 1 ,'CDEFGHJKLMNPQRSTUVWXX'[int((coordinates[1] + 80) / 8)]


def get_epsg(row,epsg):

    #search = ciop.search(end_point=reference, params=[], output_fields='self,enclosure,wkt')
    #x=loads(search[0]['wkt']).centroid.coords
    
    epsg_codes = dict()
    
    if epsg =='None':
    
        x=loads(row['wkt']).centroid.coords
        coord=x[0]
        z, l = zone(coord)
        # The EPSG code is 32600+zone for positive latitudes and 32700+zone for negatives.
        if coord[0]>=0 :
            EPSG='EPSG:326'+str(z)
        else:
            EPSG='EPSG:327'+str(z)
            
        epsg_codes['epsg'] = EPSG
        
    else:
        epsg_codes['epsg'] = epsg
    
    return pd.Series(epsg_codes)


def pre_process(gpt_path, products, aoi, resolution='10.0', polarization=None, orbit_type=None, dem_type='SRTM 1Sec HGT',show_graph=False):

    #mygraph = GraphProcessor()
    
    for index, product in products.iterrows():

        #aoi_subset = get_inteserction_aoi_prod(aoi,products['wkt'][index])
        mygraph = GraphProcessor(gpt_path)
        
        operator = 'Read'
        parameters = get_operator_default_parameters(operator)
        node_id = 'Read-{0}'.format(index)
        source_node_id = ''
        parameters['file'] = product.local_path 
        mygraph.add_node(node_id,
                         operator, 
                         parameters,
                         source_node_id)

        source_node_id = node_id
        
        operator = 'Apply-Orbit-File'

        parameters = get_operator_default_parameters(operator)
        
        if orbit_type == 'Restituted':
        
            parameters['orbitType'] = 'Sentinel Restituted (Auto Download)'
            

        node_id = 'Apply-Orbit-File-{0}'.format(index)
        mygraph.add_node(node_id, 
                         operator, 
                         parameters, 
                         source_node_id)

        source_node_id = node_id

        operator = 'ThermalNoiseRemoval'
        node_id = 'ThermalNoiseRemoval-{0}'.format(index)
        parameters = get_operator_default_parameters(operator)
        
        if polarization is not None:
            parameters['selectedPolarisations'] = polarization
            
        mygraph.add_node(node_id,
                         operator,
                         parameters,
                         source_node_id)

        source_node_id = node_id

        operator = 'Calibration'
        node_id = 'Calibration-{0}'.format(index)
        parameters = get_operator_default_parameters(operator)

        parameters['auxFile'] = 'Product Auxiliary File'
        parameters['outputImageInComplex'] = 'false'
        parameters['outputImageScaleInDb'] = 'false'
        parameters['createGammaBand'] = 'false'
        parameters['createBetaBand'] = 'false'
        parameters['selectedPolarisations'] = ''
        parameters['outputSigmaBand'] = 'true'
        parameters['outputGammaBand'] = 'false'
        parameters['outputBetaBand'] = 'false'
        
        if polarization is not None:
            
            parameters['selectedPolarisations'] = polarization
        
        mygraph.add_node(node_id,
                         operator,
                         parameters,
                         source_node_id)

        source_node_id = node_id
        
        ####Burst Merging 
        operator = 'TOPSAR-Deburst'
        parameters = get_operator_default_parameters(operator)
        node_id = 'TOPSAR-Deburst'
        
        mygraph.add_node(node_id,
                         operator,
                         parameters,
                         source_node_id)

        source_node_id = node_id
        

        operator = 'Terrain-Correction'
    
        node_id = 'Terrain-Correction-{0}'.format(index)
    
        parameters = get_operator_default_parameters(operator)

        parameters['mapProjection'] = product.epsg #'AUTO:42001'
        parameters['pixelSpacingInMeter'] = resolution           
        parameters['nodataValueAtSea'] = 'true'
        #parameters['demName'] = 'SRTM 1Sec HGT'
        #parameters['demName'] = 'SRTM 3Sec'
        parameters['demName'] = dem_type
        
        mygraph.add_node(node_id,
                         operator,
                         parameters,
                         source_node_id)
        
        source_node_id = node_id
        
        operator = 'Subset'
        
        node_id = 'Subset-{0}'.format(index)
        
        parameters = get_operator_default_parameters(operator)
        parameters['geoRegion'] = aoi
        parameters['copyMetadata'] = 'true'


        mygraph.add_node(node_id,
                         operator,
                         parameters,
                         source_node_id)

        source_node_id = node_id

        
        operator = 'Write'

        parameters = get_operator_default_parameters(operator)

        parameters['file'] = product.identifier 
        parameters['formatName'] = 'BEAM-DIMAP'

        node_id = 'Write-{0}'.format(index)

        mygraph.add_node(node_id,
                         operator,
                         parameters,
                         source_node_id)
       
        if show_graph: 
            mygraph.view_graph()
        
        mygraph.run()



def speckle_filter(gpt_path, products, working_dir='.', show_graph=True):
    
    mygraph = GraphProcessor(gpt_path)
    
    operator = 'Read'
    parameters = get_operator_default_parameters(operator)
    node_id_0 = 'Read-0'
    source_node_id = ''
    parameters['file'] = os.path.join(working_dir,'{}.dim'.format(products.identifier.values[0]))
    
    mygraph.add_node(node_id_0,
                     operator, 
                     parameters,
                     source_node_id)

    
    operator = 'Read'
    parameters = get_operator_default_parameters(operator)
    node_id_1 = 'Read-1'
    source_node_id = ''
    parameters['file'] = os.path.join(working_dir,'{}.dim'.format(products.identifier.values[1]))
    
    mygraph.add_node(node_id_1,
                     operator, 
                     parameters,
                     source_node_id)
    
    
    
    sources = dict()
    sources['master'] = node_id_0
    sources['slave'] = node_id_1
    
    operator = 'Collocate'
      
    parameters = get_operator_default_parameters(operator)
    
    node_id = 'Collocate'

    parameters['targetProductName'] = '_collocated'
    parameters['targetProductType'] = 'COLLOCATED'
    parameters['renameMasterComponents'] = 'true'
    parameters['renameSlaveComponents'] = 'true'
    parameters['masterComponentPattern']='before'
    parameters['slaveComponentPattern'] = 'after'
    parameters['resamplingType'] = 'NEAREST_NEIGHBOUR'
    
    
    mygraph.add_node(node_id,
                     operator, 
                     parameters,
                     sources)

    source_node_id = node_id
    
    operator = 'Multi-Temporal-Speckle-Filter'

    parameters = get_operator_default_parameters(operator)
    
    node_id = 'Multi-Temporal-Speckle-Filter'

    parameters['sourceBands'] = 'before,after'
    parameters['filter'] = 'Lee Sigma'
    parameters['filterSizeX'] = '3'
    parameters['filterSizeY'] = '3'
    parameters['dampingFactor'] = '2'
    parameters['estimateENL'] = 'true'
    parameters['enl'] = '1.0'
    parameters['numLooksStr'] = '1'
    parameters['windowSize'] = '7x7'
    parameters['targetWindowSizeStr'] = '3x3'
    parameters['sigmaStr'] = '0.9'
    parameters['anSize'] = '50'
    
    mygraph.add_node(node_id,
                     operator, 
                     parameters,
                     source_node_id)
    
    
    source_node_id = node_id
    
    operator = 'Write'

    parameters = get_operator_default_parameters(operator)

    parameters['file'] = 'mtsf'
    parameters['formatName'] = 'BEAM-DIMAP'

    node_id = 'Write'

    mygraph.add_node(node_id,
                     operator,
                     parameters,
                     source_node_id)

    if show_graph: 
        mygraph.view_graph()

    mygraph.run()


def create_stack(gpt_path, products, working_dir='.', show_graph=True):
    
    mygraph = GraphProcessor(gpt_path)
    
    operator = 'ProductSet-Reader'
    parameters = get_operator_default_parameters(operator)
    
    parameters['fileList'] = ','.join([ '{}.dim'.format(n) for n in products.identifier.values])
    
    node_id = 'ProductSet-Reader'
    source_node_id = ''
    
    #parameters['file'] = product.local_path 
    mygraph.add_node(node_id,
                     operator, 
                     parameters,
                     '')

    source_node_id = node_id
    
    operator = 'CreateStack'
    parameters = get_operator_default_parameters(operator)
    node_id = 'CreateStack'
        
    parameters['extent'] = 'Minimum'
    parameters['resamplingType'] = 'BICUBIC_INTERPOLATION'
    
    mygraph.add_node(node_id,
                     operator, 
                     parameters,
                     source_node_id)

    source_node_id = node_id
    
    operator = 'Write'

    parameters = get_operator_default_parameters(operator)

    parameters['file'] = 'stack'
    parameters['formatName'] = 'BEAM-DIMAP'

    node_id = 'Write'

    mygraph.add_node(node_id,
                     operator,
                     parameters,
                     source_node_id)

    if show_graph: 
        mygraph.view_graph()

    mygraph.run()
    
def list_bands(product, working_dir='.'):
    
    reader = ProductIO.getProductReader('BEAM-DIMAP')
    product = reader.readProductNodes(os.path.join(working_dir,product), None)

    return list(product.getBandNames())


def change_detection(gpt_path, input_product, output_product, expression, show_graph=False):
    
    mygraph = GraphProcessor(gpt_path)
    
    operator = 'Read'

    parameters = get_operator_default_parameters(operator)

    node_id = 'Read'

    source_node_id = ''

    parameters['file'] = input_product

    mygraph.add_node(node_id, operator, parameters, source_node_id)

    source_node_id = node_id 
    
    operator = 'BandMaths'

    parameters = get_operator_default_parameters(operator)

    bands = '''<targetBands>
        <targetBand>
          <name>change_detection</name>
          <type>float32</type>
          <expression>{}</expression>
          <description/>
          <unit/>
          <noDataValue>NaN</noDataValue>
        </targetBand>
        </targetBands>'''.format(expression)

    parameters['targetBandDescriptors'] = bands 

    node_id = 'BandMaths'

    mygraph.add_node(node_id, operator, parameters, source_node_id)

    source_node_id = node_id 
    
    operator = 'Write'

    parameters = get_operator_default_parameters(operator)

    parameters['file'] = output_product
    parameters['formatName'] = 'GeoTIFF-BigTIFF'

    node_id = 'Write'

    mygraph.add_node(node_id,
                     operator,
                     parameters,
                     source_node_id)
    
    if show_graph: 
        mygraph.view_graph()
    
    mygraph.run()
    

def convert_dim(gpt_path, input_product, show_graph=False):
    
    mygraph = GraphProcessor(gpt_path)
    
    operator = 'Read'

    parameters = get_operator_default_parameters(operator)

    node_id = 'Read'

    source_node_id = ''

    parameters['file'] = input_product

    mygraph.add_node(node_id, operator, parameters, source_node_id)

    source_node_id = node_id 
   
    operator = 'LinearToFromdB'

    node_id = 'LinearToFromdB'
    
    parameters = get_operator_default_parameters(operator)
    
    mygraph.add_node(node_id, operator, parameters, source_node_id)
    
    source_node_id = node_id 
    
    operator = 'Write'

    parameters = get_operator_default_parameters(operator)

    parameters['file'] = input_product.replace('.dim', '_db.tif')
    parameters['formatName'] = 'GeoTIFF-BigTIFF'

    node_id = 'Write'

    mygraph.add_node(node_id, operator, parameters, source_node_id)
    
    if show_graph: 
        mygraph.view_graph()
    
    mygraph.run()

    return input_product.replace('.dim', '_db.tif')
    
def cog(input_tif, output_tif):
    
    translate_options = gdal.TranslateOptions(gdal.ParseCommandLine('-co TILED=YES ' \
                                                                    '-co COPY_SRC_OVERVIEWS=YES ' \
                                                                    '-co BIGTIFF=YES ' \
                                                                    '-co COMPRESS=LZW'))

    ds = gdal.Open(input_tif, gdal.OF_READONLY)

    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
    ds.BuildOverviews('NEAREST', [2,4,8,16,32])
    
    ds = None

    ds = gdal.Open(input_tif)
    gdal.Translate(output_tif,
                   ds, 
                   options=translate_options)
    ds = None

    os.remove('{}.ovr'.format(input_tif))
    os.remove(input_tif)

def get_image_wkt(product):
    
    src = gdal.Open(product)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()

    max_x = ulx + (src.RasterXSize * xres)
    min_y = uly + (src.RasterYSize * yres)
    min_x = ulx 
    max_y = uly

    source = osr.SpatialReference()
    source.ImportFromWkt(src.GetProjection())

    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(source, target)

    result_wkt = box(transform.TransformPoint(min_x, min_y)[0],
                     transform.TransformPoint(min_x, min_y)[1],
                     transform.TransformPoint(max_x, max_y)[0],
                     transform.TransformPoint(max_x, max_y)[1]).wkt
    
    return result_wkt


def create_composite(input_products, output_product, band_expressions):
    
    BandMathX = otbApplication.Registry.CreateApplication("BandMathX")

    BandMathX.SetParameterStringList('il', input_products)
    BandMathX.SetParameterString('out', 'temp_red_green_blue.tif')
    BandMathX.SetParameterString('exp', ';'.join(band_expressions))

    BandMathX.ExecuteAndWriteOutput()

    Convert = otbApplication.Registry.CreateApplication('DynamicConvert')

    Convert.SetParameterString('in', 'temp_red_green_blue.tif')
    Convert.SetParameterString('out', output_product)
    Convert.SetParameterString('type', 'linear')
    Convert.SetParameterString('channels', 'rgb')

    Convert.ExecuteAndWriteOutput()

    os.remove('temp_red_green_blue.tif')
    
    return output_product

def create_mask(in_composite, out_mask):
    
    #gdal_calc.py --calc="logical_and(logical_and(A==255, B==0), C==0)" -A $1 --A_band=1 -B $1 --B_band=2 -C $1 --C_band=3 --outfile=${1::-8}.mask.tif
    #command = '/opt/anaconda/envs/env_ewf_satcen_03_03_01/bin/gdal_calc.py --calc="logical_and(logical_and(A==255, B==0), C==0)" -A {0} --A_band=1 -B {0} --B_band=2 -C {0} --C_band=3 --outfile={1}'.format(in_composite, out_mask)

    # Run the command. os.system() returns value zero if the command was executed succesfully
    #out = os.system(command)
    #if out !=0:
    #    print('ERROR gdla_calc out: {}'.format(out))
    
    calc_exp="logical_and(logical_and(A==255, B==0), C==0)"
    
    gdalCalc(calc=calc_exp, A=in_composite, A_band=1, B=in_composite, B_band=2, C=in_composite, C_band=3, outfile=out_mask)
    
    
def create_rbb(in_rgb, out_rbb):
    
    #gdal_translate -ot UInt16 -a_nodata 256 ${1::-14}RED-BLUE.rgb.tif ${1::-8}.acd.tif -co COMPRESS=LZW -b 1 -b 3 -b 3
    
    translate_options = gdal.TranslateOptions(gdal.ParseCommandLine('-co COMPRESS=LZW '\
                                                                    '-ot UInt16 ' \
                                                                    '-a_nodata 256 ' \
                                                                    '-b 1 -b 3 -b 3 '))
                                                                    
    ds = gdal.Open(in_rgb, gdal.OF_READONLY)

    gdal.Translate(out_rbb, 
                   ds, 
                   options=translate_options)
    
    
def write_tif(layer, output_tif, width, height, input_geotransform, input_georef):

    driver = gdal.GetDriverByName('GTiff')

    ds = driver.Create(output_tif, 
                       width, 
                       height, 
                       1, 
                       gdal.GDT_UInt16) 
        
    ds.SetGeoTransform(input_geotransform)
    ds.SetProjection(input_georef)
    ds.GetRasterBand(1).WriteArray(layer)
        
    ds.FlushCache()
    ds = None

    
def create_mask_2(input_tif, output_tif):
    
    ds = gdal.Open(input_tif, gdal.OF_READONLY)
    input_geotransform = ds.GetGeoTransform()
    input_georef = ds.GetProjectionRef()
    
    w = ds.GetRasterBand(1).XSize
    h = ds.GetRasterBand(1).YSize
    band1 = ds.GetRasterBand(1).ReadAsArray(0,0,w,h)
    band2 = ds.GetRasterBand(2).ReadAsArray(0,0,w,h)
    band3 = ds.GetRasterBand(3).ReadAsArray(0,0,w,h)

    out = np.zeros((h, w))
    out[(np.where((band1==255) & (band2==0) & band3==0))] = 1
    
    ds = None
    
    write_tif(out, output_tif, w, h, input_geotransform, input_georef)