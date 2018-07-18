/////////////////////////////////////////////////////////////////////
//PeatLandProbabilityInputs - GEE script                            /
//Written by:                                                       /
////Evan R. DeLancey                                                /
////GIS Land-use analyst Alberta Biodiversity Monitoring Institute  /
////2017-08-24                                                      /
////edelance@ualberta.ca                                            /
/////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
//Define Alberta shape and projection file
///////////////////////////////////////////////////////////////////////////////
var PrjFile = ee.Image('users/abmigc/ProcessingUnits');
var prj = PrjFile.projection();
var AB = ee.FeatureCollection('users/abmigc/Alberta_prjGEE');
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////








///////////////////////////////////////////////////////////////////////////////
//Define Sentinel-1 functions
///////////////////////////////////////////////////////////////////////////////

//Angle correction and add NDPOL
function S1Prep(image) {
  var first = image.addBands(image.select('VV').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VV_Gamma'));
  var second = first.addBands(image.select('VH').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VH_Gamma'));
  return second.addBands(second.normalizedDifference(['VH_Gamma', 'VV_Gamma']).rename('NDPOL'));
}

//Sigma Lee filter
function toNatural(img) {
  return ee.Image(10.0).pow(img.select(0).divide(10.0));
}
function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}
// The RL speckle filter from https://code.earthengine.google.com/2ef38463ebaf5ae133a478f173fd0ab5
// by Guido Lemoine
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}

//angle masking
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.53993));
};
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.53993));
};
var maskAng3 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var maskAng4 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(44.73993));
};
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////

 

///////////////////////////////////////////////////////////////////
// START Get Sentinel-1 image stack for 2016 and 2017
///////////////////////////////////////////////////////////////////
var s1_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2016-07-01', '2016-07-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
var s1_1 = s1_1.map(maskAng1);
var s1_1 = s1_1.map(maskAng2);

var s1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2016-08-01', '2016-08-10')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
var s1_2 = s1_2.map(maskAng1);
var s1_2 = s1_2.map(maskAng2);

var s1_2b = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2016-08-11', '2016-08-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
var s1_2b = s1_2b.map(maskAng1);
var s1_2b = s1_2b.map(maskAng2);

var s1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2017-05-15', '2017-08-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
var s1_3 = s1_3.map(maskAng3);
var s1_3 = s1_3.map(maskAng4);

var s1 = ee.ImageCollection(s1_1.merge(s1_2));
var s1 = ee.ImageCollection(s1.merge(s1_3));
var s1 = ee.ImageCollection(s1.merge(s1_2b));
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////
//Map functions over Sentinel-1 image stack
///////////////////////////////////////////////////////////////////////////////
var s1 = s1.map(S1Prep);
var VH = s1.select(['VH_Gamma']);
var NDPOL = s1.select(['NDPOL']);
var VH = VH.map(toNatural);
var VH = VH.map(RefinedLee);
var VH = VH.map(toDB);
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////
//Spatialy filter NDPOL
///////////////////////////////////////////////////////////////////////////////
var boxcar = ee.Kernel.circle({
  radius: 3, units: 'pixels', normalize: true
});
function fltr(image) {
  return image.convolve(boxcar);
}
var NDPOL = NDPOL.map(fltr);
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////
//Take temporal mean of picel stack and add layers to map
///////////////////////////////////////////////////////////////////////////////
var VH = VH.mean();
var NDPOL = NDPOL.mean();
Map.addLayer(VH, {min:-28, max:-8}, 'VH');
Map.addLayer(NDPOL, {min:0, max:0.5}, 'NPOL');
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////








///////////////////////////////////////////////////////////////////////////////
//Reproject and clip Sentinel-1 variables to Alberta
///////////////////////////////////////////////////////////////////////////////
var VH = VH.reproject(prj, null, 10).clip(AB);
var NDPOL = NDPOL.reproject(prj, null, 10).clip(AB);
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////
//Export Sentinel-1 variables to drive
/////////////////////////////////////////////////////////////////
Export.image.toDrive({
  image: VH,
  description: 'VH',
  scale: 10,
  region: AB,
  maxPixels: 3E10
});
Export.image.toDrive({
  image: NDPOL,
  description: 'NDPOL',
  scale: 10,
  region: AB,
  maxPixels: 3E10
});
///////////////////////////////////////////////////////////////
// END
///////////////////////////////////////////////////////////////












///////////////////////////////////////////////////////////////
// Sentinel-2 functions
///////////////////////////////////////////////////////////////

//mask clouds
function maskCloud(image) {
  var QA60 = image.select(['QA60']);
  var B1 = image.select(['B1']).gt(1500);
  var mask1 = image.updateMask(QA60.lt(1));
  return mask1.updateMask(B1.lt(1));
}

//add indices
function addIndices(image) {
 var a = image.addBands(image.normalizedDifference(['B8', 'B4']).rename('NDVI'));
 var b = a.addBands(a.normalizedDifference(['B3', 'B8']).rename('NDWI'));
 var c = b.addBands(b.expression(
      '(B8 / B2) - (B8 / B3)', {
        'B8': image.select(['B8']),
        'B2': image.select(['B2']),
        'B3': image.select(['B3'])
      }
    ).rename('ARI'));
  var d = c.addBands(c.expression(
      '(B4 - B2) / B5', {
        'B4': image.select(['B4']),
        'B2': image.select(['B2']),
        'B5': image.select(['B5'])
      }
    ).rename('PSRI'));
  return d.addBands(d.expression(
      '705 + 35*((((RED + RE3)/2) - RE1) / (RE2 - RE1))', {
        'RE1': image.select(['B5']),
        'RE2': image.select(['B6']),
        'RE3': image.select(['B7']),
        'RED' : image.select(['B4'])
      }
    ).rename('REIP'));
}

///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////









///////////////////////////////////////////////////////////////
//Get Sentinel-1 image stack
///////////////////////////////////////////////////////////////
var S2 = ee.ImageCollection('COPERNICUS/S2')
  .filterDate('2016-05-15', '2016-08-31')
  .filterBounds(AB);
var S2_1 = ee.ImageCollection('COPERNICUS/S2')
  .filterDate('2017-05-15', '2017-08-31')
  .filterBounds(AB);
var S2 = ee.ImageCollection(S2.merge(S2_1));
///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////









///////////////////////////////////////////////////////////////
//Map functions over image stack
///////////////////////////////////////////////////////////////
var S2 = S2.map(maskCloud);
var S2 = S2.map(addIndices);
///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////









///////////////////////////////////////////////////////////////
//Get median value of each index and band and add to map
///////////////////////////////////////////////////////////////
var B2 = S2.select(['B2']).median();
var B3 = S2.select(['B3']).median();
var B4 = S2.select(['B4']).median();
var B8 = S2.select(['B8']).median();
var NDVI = S2.select(['NDVI']).median();
var NDWI = S2.select(['NDWI']).median();
var ARI = S2.select(['ARI']).median();
var PSRI = S2.select(['PSRI']).median();
var REIP = S2.select(['REIP']).median();

var colorbrewer = require('users/gena/packages:colorbrewer');

Map.addLayer(NDVI, {min:-0.5, max:0.9, palette: colorbrewer.Palettes.RdYlGn[11]}, 'NDVI');
Map.addLayer(NDWI, {min:-1, max:1, palette: colorbrewer.Palettes.Blues[9]}, 'NDWI');
Map.addLayer(ARI, {min:-1, max:0.3, palette: colorbrewer.Palettes.PRGn[11]}, 'ARI');
Map.addLayer(PSRI, {min:-1, max:1, palette: colorbrewer.Palettes.PRGn[11]}, 'PSRI');
Map.addLayer(REIP, {min:715, max:730, palette: colorbrewer.Palettes.PRGn[11]}, 'REIP');
///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////








///////////////////////////////////////////////////////////////
//Reproject indices and clip to Alberta
///////////////////////////////////////////////////////////////
var B2 = B2.reproject(prj, null, 10).clip(AB);
var B3 = B3.reproject(prj, null, 10).clip(AB);
var B4 = B4.reproject(prj, null, 10).clip(AB);
var B8 = B8.reproject(prj, null, 10).clip(AB);
var NDVI = NDVI.reproject(prj, null, 10).clip(AB);
var NDWI = NDWI.reproject(prj, null, 10).clip(AB);
var ARI = ARI.reproject(prj, null, 10).clip(AB);
var PSRI = PSRI.reproject(prj, null, 10).clip(AB);
var REIP = REIP.reproject(prj, null, 10).clip(AB);
///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////










///////////////////////////////////////////////////////////////
//Export Sentinel-2 indices to drive
///////////////////////////////////////////////////////////////
Export.image.toDrive({
  image: B2,
  description: 'B2',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: B3,
  description: 'B3',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: B4,
  description: 'B4',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: B8,
  description: 'B8',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: NDWI,
  description: 'NDWI',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: NDVI,
  description: 'NDVI',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: ARI,
  description: 'ARI',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: PSRI,
  description: 'PSRI',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
Export.image.toDrive({
  image: REIP,
  description: 'REIP',
  scale: 10,
  region: AB,
  maxPixels: 10E10
});
///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////
//Get SRTM DEM, smooth, and export to drive
///////////////////////////////////////////////////////////////
var srtm = ee.Image('USGS/SRTMGL1_003');
var srtm = srtm.reproject(prj, null, 10).clip(AB);
var boxcar = ee.Kernel.circle({
  radius: 7, units: 'pixels', normalize: true
});
var srtm = srtm.convolve(boxcar);

Export.image.toDrive({
  image: srtm,
  description: 'SRTM',
  scale: 10,
  region: AB,
  maxPixels: 9E10
});
///////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////
//END script
///////////////////////////////////////////////////////////////
