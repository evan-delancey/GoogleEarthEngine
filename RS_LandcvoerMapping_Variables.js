/////////////////////////////////////////////////////////////////////
//Boreal Probability of Wet Area Get Variables                      /
//Written by:                                                       /
////Evan R. DeLancey                                                /
////GIS Land-use analyst Alberta Biodiversity Monitoring Institute  /
////2017-08-24                                                      /
////edelance@ualberta.ca                                            /
/////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////
// START Get Sentinel-1 VH and NPOL
///////////////////////////////////////////////////////////////////
//get sentinel-1 images filtered by date, polarization, resolution, and orbit
var s1_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2016-07-01', '2016-07-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);

//mask out edges of images by using angle
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.53993));
};
var s1_1 = s1_1.map(maskAng1);

//mask out edges of images by using angle
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.53993));
};
var s1_1 = s1_1.map(maskAng2);

var ClipImgColl = function(image){
  return image.clip(July_Clip1);
};

var s1_1 = s1_1.map(ClipImgColl);

  
var s1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2016-08-01', '2016-08-10')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);

//mask out edges of images by using angle
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.53993));
};
var s1_2 = s1_2.map(maskAng1);

//mask out edges of images by using angle
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.53993));
};
var s1_2 = s1_2.map(maskAng2);

var ClipImgColl = function(image){
  return image.clip(August_Clip1);
};

var s1_2 = s1_2.map(ClipImgColl);

var s1_2b = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2016-08-11', '2016-08-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);

//mask out edges of images by using angle
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.53993));
};
var s1_2b = s1_2b.map(maskAng1);

//mask out edges of images by using angle
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.33993));
};
var s1_2b = s1_2b.map(maskAng2);

var ClipImgColl = function(image){
  return image.clip(August_Clip2);
};

var s1_2b = s1_2b.map(ClipImgColl);

var s1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2017-05-15', '2017-08-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
  
//mask out edges of images by using angle
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_3 = s1_3.map(maskAng1);

//mask out edges of images by using angle
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(44.73993));
};
var s1_3 = s1_3.map(maskAng2);


  
var s1 = ee.ImageCollection(s1_1.merge(s1_2));
var s1 = ee.ImageCollection(s1.merge(s1_3));
var s1 = ee.ImageCollection(s1.merge(s1_2b));
print(s1);



//angle correction
function toGamma0(image) {
  var vh = image.select('VH').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0));
  return vh.addBands(image.select('VV').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0)));
}



//filter out windy days
var pctWat = function(image){
  var d = image.date().format('Y-M-d');
  var wx = ee.ImageCollection('NOAA/CFSV2/FOR6H')
    .filterDate(d);
  var vWind = wx.select(['v-component_of_wind_height_above_ground']);
  var a = vWind.max();
  var uWind = wx.select(['u-component_of_wind_height_above_ground']);
  var b = uWind.max();
  var a = a.pow(2);
  var b = b.pow(2);
  var ab = a.add(b);
  var ws = ab.sqrt();
  var ws = ws.multiply(3.6);
  return image.updateMask(ws.lt(12));
};
var s1 = s1.map(pctWat);



//add in difference between vv and vh
var adddiff = function(image) {
  return image.addBands(image.expression(
      '(VH - VV) / (VH + VV)', {
        'VH': image.select(['VH']),
        'VV': image.select(['VV'])
      }
    ));
};

var s1 = s1.map(toGamma0);

var s1 = s1.map(adddiff);
var VH = s1.select(['VH']);
var NPOL = s1.select(['VH_1']);

////////////////////////////////////////////////////
//Sigma Lee filter
////////////////////////////////////////////////////
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

var VH = VH.map(toNatural);
var VH = VH.map(RefinedLee);
var VH = VH.map(toDB);
/////////////////////////////////////////////////

var VH = VH.mean();

var boxcar = ee.Kernel.circle({
  radius: 3, units: 'pixels', normalize: true
});

var fltr = function(image) {
  return image.convolve(boxcar);
};

var NPOL = NPOL.map(fltr);
var NPOL = NPOL.mean();


//Map.addLayer(VH, {min:-28, max:-8}, 'VH');
//Map.addLayer(NPOL, {min:0, max:0.5}, 'NPOL');
///////////////////////////////////////////////////////////////////
// END Get Sentinel-1 VH and NPOL
///////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////
// START Get Sentinel-1 VVSD
///////////////////////////////////////////////////////////////////
var s1_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2016-05-15', '2016-07-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
          return img.clip(img.geometry().buffer(-12000));
  });

var s1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2016-08-01', '2016-08-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
          return img.clip(img.geometry().buffer(-14000));
  });
  
var s1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2017-05-15', '2017-08-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
          return img.clip(img.geometry().buffer(-1000));
  });
  
//mask out edges of images by using angle
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_3 = s1_3.map(maskAng1);

//mask out edges of images by using angle
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.23993));
};
var s1_3 = s1_3.map(maskAng2);

var s1_4 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filterDate('2015-05-15', '2015-08-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
          return img.clip(img.geometry().buffer(-10000));
  });
  
var S1 = ee.ImageCollection(s1_1.merge(s1_2));
var S1 = ee.ImageCollection(S1.merge(s1_3));
var S1 = ee.ImageCollection(S1.merge(s1_4));

function toGamma01(image) {
  return image.select('VV').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0));
}


  
var S1 = S1.map(toGamma01);
var S1 = S1.map(toNatural);
var S1 = S1.map(RefinedLee);
var S1 = S1.map(toDB);
print(S1);
var VV = S1.select(['sum']);



var VVSD = VV.reduce(ee.Reducer.stdDev());


Map.addLayer(VVSD, {min:1, max:5}, 'VVSD');

///////////////////////////////////////////////////////////////////
// END Get Sentinel-1 VVSDL
///////////////////////////////////////////////////////////////////


var prj = image.projection();
var VH = VH.resample('bicubic').reproject(prj, null, 10);
var VH = VH.clip(geometry);
var NPOL = NPOL.resample('bicubic').reproject(prj, null, 10);
var NPOL = NPOL.clip(geometry);
var VVSD = VVSD.resample('bicubic').reproject(prj, null, 10);
var VVSD1 = VVSD.clip(G1);
var VVSD2 = VVSD.clip(G2);
var VVSD3 = VVSD.clip(G3);

//////////////////////////////////////////////////////////////////
// START export Sentinel-1 variables
/////////////////////////////////////////////////////////////////

Export.image.toDrive({
  image: VH,
  description: 'VH',
  scale: 10,
  region: geometry,
  maxPixels: 3E10
  
});

Export.image.toDrive({
  image: NPOL,
  description: 'NPOL',
  scale: 10,
  region: geometry,
  maxPixels: 3E10
  
});

Export.image.toDrive({
  image: VVSD1,
  description: 'VVSD_1',
  scale: 10,
  region: Ge1,
  maxPixels: 3E10
  
});

Export.image.toDrive({
  image: VVSD2,
  description: 'VVSD_2',
  scale: 10,
  region: Ge2,
  maxPixels: 3E10
  
});

Export.image.toDrive({
  image: VVSD3,
  description: 'VVSD_3',
  scale: 10,
  region: Ge3,
  maxPixels: 3E10
  
});

Map.addLayer(G3, {}, 'G1');


///////////////////////////////////////////////////////////////
// END export Sentinel-1 varibles
///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
// START Get Sentinel-2 and cloud removal
///////////////////////////////////////////////////////////////

var S2 = ee.ImageCollection('COPERNICUS/S2')
  //filter start and end date
  .filterDate('2016-05-15', '2016-08-31')
  //filter according to drawn boundary
  .filterBounds(geometry);
  
//Get Sentinel-2 data
var S2_1 = ee.ImageCollection('COPERNICUS/S2')
  //filter start and end date
  .filterDate('2017-05-15', '2017-08-31')
  //filter according to drawn boundary
  .filterBounds(geometry);

var S2 = ee.ImageCollection(S2.merge(S2_1));
print(S2);

//mask cloud from built of cloud mask
var maskcloud1 = function(image) {
  var QA60 = image.select(['QA60']);
  return image.updateMask(QA60.lt(1));
};
var S2 = S2.map(maskcloud1);

//mask further cloud with B1 threshold
var maskcloud2 = function(image) {
  var B1 = image.select(['B1']);
  var B11 = image.select(['B11']);
  var B12 = image.select(['B12']);
  var bin = B1.gt(1500);
  return image.updateMask(bin.lt(1));
};
var S2 = S2.map(maskcloud2);


// add NDVI to bands
var addNDVI = function(image) {
  return image.addBands(image.normalizedDifference(['B8', 'B4']));
};
var S2 = S2.map(addNDVI);
var NDVI = S2.select(['nd']);
var NDVI = NDVI.median();

// add NDWI to bands
var addNDWI = function(image) {
  return image.addBands(image.normalizedDifference(['B3', 'B8']));
};
var S2 = S2.map(addNDWI);
var NDWI = S2.select(['nd_1']);
var NDWI = NDWI.median();


// add ARI to bands
var addARI = function(image) {
  return image.addBands(image.expression(
      '(B8 / B2) - (B8 / B3)', {
        'B8': image.select(['B8']),
        'B2': image.select(['B2']),
        'B3': image.select(['B3'])
      }
    ));
};
var S2 = S2.map(addARI);
var ARI = S2.select(['B8_1']);
var ARI = ARI.median();

var addPSRI = function(image) {
  return image.addBands(image.expression(
      '(B4 - B2) / B5', {
        'B4': image.select(['B4']),
        'B2': image.select(['B2']),
        'B5': image.select(['B5'])
      }
    ));
};
var S2 = S2.map(addPSRI);
var PSRI = S2.select(['B4_1']);
var PSRI = PSRI.median();

var addREIP = function(image) {
  return image.addBands(image.expression(
      '705 + 35*((((RED + RE3)/2) - RE1) / (RE2 - RE1))', {
        'RE1': image.select(['B5']),
        'RE2': image.select(['B6']),
        'RE3': image.select(['B7']),
        'RED' : image.select(['B4'])
      }
    ));
};
var S2 = S2.map(addREIP);
print(S2);
var REIP = S2.select(['constant']);
var REIP = REIP.median();

var ndvi_pal = ['#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#a6d96a', '#66bd63', '#1a9850'];
Map.addLayer(NDVI, {min:-0.5, max:0.9, palette: ndvi_pal}, 'NDVI');
var ndwi_pal = ['#ece7f2', '#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0', '#045a8d', '#023858'];
Map.addLayer(NDWI, {min:-1, max:1, palette: ndwi_pal}, 'NDWI');
var ari_pal = ['#1b7837', '#5aae61', '#a6dba0', '#d9f0d3', '#e7d4e8', '#c2a5cf', '#9970ab', '#762a83'];
Map.addLayer(ARI, {min:-1, max:0.3, palette: ari_pal}, 'ARI');
Map.addLayer(PSRI, {min:-1, max:1, palette:ari_pal}, 'PSRI');
Map.addLayer(REIP, {min:715, max:730, palette:ari_pal}, 'REIP');

var B2 = S2.select(['B2']);
var B2 = B2.median();

var B3 = S2.select(['B3']);
var B3 = B3.median();

var B4 = S2.select(['B4']);
var B4 = B4.median();

var B8 = S2.select(['B8']);
var B8 = B8.median();

var NDVI = NDVI.reproject(prj, null, 10);
var NDVI = NDVI.clip(geometry);

var NDWI = NDWI.reproject(prj, null, 10);
var NDWI = NDWI.clip(geometry);

var ARI = ARI.reproject(prj, null, 10);
var ARI = ARI.clip(geometry);

var PSRI = PSRI.reproject(prj, null, 10);
var PSRI = PSRI.clip(geometry);

var REIP = REIP.reproject(prj, null, 10);
var REIP = REIP.clip(geometry);

var B2 = B2.reproject(prj, null, 10);
var B2 = B2.clip(geometry);

var B3 = B3.reproject(prj, null, 10);
var B3 = B3.clip(geometry);

var B4 = B4.reproject(prj, null, 10);
var B4 = B4.clip(geometry);

var B8 = B8.reproject(prj, null, 10);
var B8 = B8.clip(geometry);


Export.image.toDrive({
  image: REIP,
  description: 'REIP',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: PSRI,
  description: 'PSRI',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: ARI,
  description: 'ARI',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: NDWI,
  description: 'NDWI',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: NDVI,
  description: 'NDVI',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: B2,
  description: 'B2',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: B3,
  description: 'B3',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: B4,
  description: 'B4',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: B8,
  description: 'B8',
  scale: 10,
  region: geometry,
  maxPixels: 10E10
  
});
