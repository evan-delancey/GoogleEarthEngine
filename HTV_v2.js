/////////////////////////////////////////////////////////////////////////////////////////
//Hydro Temporal Varriation version 2
//written by:
///Evan R. DeLancey
///GIS Land-use Analyst Alberta Biodiversity Monitoring Insitute
////////////////////////////////////////////////////////////////////////////////////////


//#####################################################################################
//#####################################################################################
///////////////////////////////////////////////////////////////////////////////////////
//1.Run Grassland HTV Algorithm
//1.START
///////////////////////////////////////////////////////////////////////////////////////
//#####################################################################################
//#####################################################################################


//Define grasslands zone
var grasslands = ee.FeatureCollection('users/abmigc/Grasslands_Zone1');

///////////////////////////////////////////////////////////////////////////////////////
//1.1 Get Sentinel-1 image stack
//1.1 START
///////////////////////////////////////////////////////////////////////////////////////

//2014 data
var s1_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(grasslands)
  .filterDate('2014-04-01', '2014-10-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-8000));
  });
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_1 = s1_1.map(maskAng1);

var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.23993));
};
var s1_1 = s1_1.map(maskAng2);



//2015 data
var s1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(grasslands)
  .filterDate('2015-04-01', '2015-10-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-1000));
  });
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_2 = s1_2.map(maskAng1);
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.23993));
};
var s1_2 = s1_2.map(maskAng2);




//2016 data
var s1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(grasslands)
  .filterDate('2016-04-01', '2016-10-31')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-10000));
  });
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_3 = s1_3.map(maskAng1);
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.23993));
};
var s1_3 = s1_3.map(maskAng2);


//get 2017 data
var s1_4 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(grasslands)
  .filterDate('2017-04-01', '2017-08-09')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-1000));
  });
  

var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_4 = s1_4.map(maskAng1);
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(44.73993));
};
var s1_4 = s1_4.map(maskAng2);


//merge 2014-2017
var s1 = ee.ImageCollection(s1_1.merge(s1_2));
var s1 = ee.ImageCollection(s1.merge(s1_3));
var s1 = ee.ImageCollection(s1.merge(s1_4));

//save image stack for later SD processing
var s1_for_SD = s1;

///////////////////////////////////////////////////////////////////////////////////////
//1.1 Get Sentinel-1 image stack
//1.1 END
///////////////////////////////////////////////////////////////////////////////////////



//####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//1.2 Filter out windy days
//1.2 START
///////////////////////////////////////////////////////////////////////////////////////


var FilterWind = function(image){
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
  return image.updateMask(ws.lt(9));
};

var s1 = s1.map(FilterWind);
var count_Grasslands = s1.select(['VV']);
var count_Grasslands = count_Grasslands.count();

///////////////////////////////////////////////////////////////////////////////////////
//1.2 Filter out windy days
//1.2 END
///////////////////////////////////////////////////////////////////////////////////////



//####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//1.3 gamma0 angle correction
//1.3 START
///////////////////////////////////////////////////////////////////////////////////////

//angle correction
function toGamma0(image) {
  return image.select('VV').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0));
}
var s1 = s1.map(toGamma0);

///////////////////////////////////////////////////////////////////////////////////////
//1.3 gamma0 angle correction
//1.3 START
///////////////////////////////////////////////////////////////////////////////////////



//####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//1.4 Sigma Lee filter
//1.4 START
///////////////////////////////////////////////////////////////////////////////////////

var s1 = s1.select(['VV']);
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

var s1 = s1.map(toNatural);
var s1 = s1.map(RefinedLee);
var s1 = s1.map(toDB);

///////////////////////////////////////////////////////////////////////////////////////
//1.4 Sigma Lee filter
//1.4 END
///////////////////////////////////////////////////////////////////////////////////////



//####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//1.5 Calculate HTV from pixel stack
//1.5 START
///////////////////////////////////////////////////////////////////////////////////////

var pctWat = function(image){
  return image.lt(-17.5);
};

var count = s1.count();
var wat = s1.map(pctWat);
var wat = wat.sum();

var PercentWater = wat.divide(count).multiply(100);
var PercentWater = PercentWater.toInt();
//Map.addLayer(PercentWater, {min:0, max:100}, 'watpct');

var prj = image.projection();
var PercentWater = PercentWater.reproject(prj, null, 10);
var PercentWater = PercentWater.clip(grasslands);


///////////////////////////////////////////////////////////////////////////////////////
//1.5 Calculate HTV from pixel stack
//1.5 END
///////////////////////////////////////////////////////////////////////////////////////



//####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//1.6 Get standard deviation of pixel stack
//1.6 START
///////////////////////////////////////////////////////////////////////////////////////

var s1_asc = ee.ImageCollection(s1_for_SD)
  .filterMetadata('orbitProperties_pass', 'equals' , 'ASCENDING');
var s1_asc = s1_asc.map(toGamma0);
var s1_asc = s1_asc.map(toNatural);
var s1_asc = s1_asc.map(RefinedLee);
var s1_asc = s1_asc.map(toDB);
var s1_asc_SD = s1_asc.reduce(ee.Reducer.stdDev());

var s1_dsc = ee.ImageCollection(s1_for_SD)
  .filterMetadata('orbitProperties_pass', 'equals' , 'DESCENDING');
var s1_dsc = s1_dsc.map(toGamma0);
var s1_dsc = s1_dsc.map(toNatural);
var s1_dsc = s1_dsc.map(RefinedLee);
var s1_dsc = s1_dsc.map(toDB);
var s1_dsc_SD = s1_dsc.reduce(ee.Reducer.stdDev());

var s1_SD = s1_asc_SD.add(s1_dsc_SD).divide(2);
var s1_SD = s1_SD.reproject(prj, null, 10);
var s1_SD = s1_SD.clip(grasslands);

var s1_mean =  ee.ImageCollection(s1_asc.merge(s1_dsc));
var s1_mean = s1_mean.mean();
var s1_mean = s1_mean.reproject(prj, null, 10);
var s1_mean = s1_mean.clip(grasslands);

var count_Grasslands_50 = count_Grasslands.reproject(prj, null, 50);
var count_Grasslands_50 = count_Grasslands_50.clip(grasslands);

var count_Grasslands_100 = count_Grasslands.reproject(prj, null, 100);
var count_Grasslands_100 = count_Grasslands_100.clip(grasslands);





///////////////////////////////////////////////////////////////////////////////////////
//1.6 Get standard deviation of pixel stack
//1.6 END
///////////////////////////////////////////////////////////////////////////////////////


//####################################################################################


///////////////////////////////////////////////////////////////////////////////////////
//1.7 Export Images to Google Drive
//1.7 START
///////////////////////////////////////////////////////////////////////////////////////

Export.image.toDrive({
  image: s1_SD,
  description: 'S1_SD_Grasslands',
  scale: 10,
  region: grasslands,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: PercentWater,
  description: 'HTV_Grasslands',
  scale: 10,
  region: grasslands,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: s1_mean,
  description: 's1_mean',
  scale: 10,
  region: grasslands,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: count_Grasslands_50,
  description: 'count_Grasslands_50',
  scale: 50,
  region: grasslands,
  maxPixels: 10E10
  
});
Export.image.toDrive({
  image: count_Grasslands_100,
  description: 'count_Grasslands_100',
  scale: 100,
  region: grasslands,
  maxPixels: 10E10
  
});

///////////////////////////////////////////////////////////////////////////////////////
//1.7 Export Images to Google Drive
//1.7 START
///////////////////////////////////////////////////////////////////////////////////////




//#####################################################################################
//#####################################################################################
///////////////////////////////////////////////////////////////////////////////////////
//1.Run Grassland HTV Algorithm
//1.END
///////////////////////////////////////////////////////////////////////////////////////
//#####################################################################################
//#####################################################################################











//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################











//#####################################################################################
//#####################################################################################
///////////////////////////////////////////////////////////////////////////////////////
//2.Run Boreal HTV Algorithm
//2.START
///////////////////////////////////////////////////////////////////////////////////////
//#####################################################################################
//#####################################################################################

//Define Boreal zone
var boreal = ee.FeatureCollection('users/abmigc/Boreal_Zone2');
var boreal_1 = ee.FeatureCollection('users/abmigc/Boreal_Zone2_1');
var boreal_2 = ee.FeatureCollection('users/abmigc/Boreal_Zone2_2');


///////////////////////////////////////////////////////////////////////////////////////
//2.1 Get Sentinel-1 image stack
//2.1 START
///////////////////////////////////////////////////////////////////////////////////////

//2015 data
var s1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(boreal)
  .filterDate('2015-05-01', '2015-09-30')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-1000));
  });
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_2 = s1_2.map(maskAng1);
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.23993));
};
var s1_2 = s1_2.map(maskAng2);


//2016 data
var s1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(boreal)
  .filterDate('2016-05-01', '2016-09-30')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-13000));
  });
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_3 = s1_3.map(maskAng1);
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(45.23993));
};
var s1_3 = s1_3.map(maskAng2);

//2017 data
var s1_4 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(boreal)
  .filterDate('2017-05-01', '2017-08-09')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filterMetadata('resolution_meters', 'equals' , 10)
  .map(function (img) {
        return img.clip(img.geometry().buffer(-1000));
  });
var maskAng1 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.gt(30.63993));
};
var s1_4 = s1_4.map(maskAng1);
var maskAng2 = function(image) {
  var ang = image.select(['angle']);
  return image.updateMask(ang.lt(44.73993));
};
var s1_4 = s1_4.map(maskAng2);


//merge 2015-2017
var s1 = ee.ImageCollection(s1_2.merge(s1_3));
var s1 = ee.ImageCollection(s1.merge(s1_4));

//save image stack for later SD processing
var s1_for_SD = s1;

///////////////////////////////////////////////////////////////////////////////////////
//2.1 Get Sentinel-1 image stack
//2.1 END
///////////////////////////////////////////////////////////////////////////////////////



//#####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//2.2 Filter out windy days
//2.2 START
///////////////////////////////////////////////////////////////////////////////////////

var FilterWind = function(image){
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
  return image.updateMask(ws.lt(9));
};
var s1 = s1.map(FilterWind);
var count_Boreal = s1.select(['VV']);
var count_Boreal = count_Boreal.count();

///////////////////////////////////////////////////////////////////////////////////////
//2.2 Filter out windy days
//2.2 END
///////////////////////////////////////////////////////////////////////////////////////



//#####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//2.3 Gamma0 angle correction
//2.3 START
///////////////////////////////////////////////////////////////////////////////////////

var s1 = s1.map(toGamma0);

///////////////////////////////////////////////////////////////////////////////////////
//2.3 Gamma0 angle correction
//2.3 END
///////////////////////////////////////////////////////////////////////////////////////



//#####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//2.4 Sigma Lee filter
//2.4 START
///////////////////////////////////////////////////////////////////////////////////////

var s1 = s1.select(['VV']);
var s1 = s1.map(toNatural);
var s1 = s1.map(RefinedLee);
var s1 = s1.map(toDB);

///////////////////////////////////////////////////////////////////////////////////////
//2.4 Sigma Lee filter
//2.4 END
///////////////////////////////////////////////////////////////////////////////////////


//#####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//2.5 Calculate HTV values
//2.5 START
///////////////////////////////////////////////////////////////////////////////////////

var pctWat = function(image){
  return image.lt(-15.1);
};

var count = s1.count();
var wat = s1.map(pctWat);
var wat = wat.sum();
var PercentWater = wat.divide(count).multiply(100);
var PercentWater = PercentWater.toInt();
var prj = image.projection();
var PercentWater = PercentWater.reproject(prj, null, 10);
var PercentWater = PercentWater.clip(boreal);


///////////////////////////////////////////////////////////////////////////////////////
//2.5 Calculate HTV values
//2.5 END
///////////////////////////////////////////////////////////////////////////////////////



//#####################################################################################



///////////////////////////////////////////////////////////////////////////////////////
//2.6 Get standard deviation of pixel stack
//2.6 START
///////////////////////////////////////////////////////////////////////////////////////

var s1_asc = ee.ImageCollection(s1_for_SD)
  .filterMetadata('orbitProperties_pass', 'equals' , 'ASCENDING');
var s1_asc = s1_asc.map(toGamma0);
var s1_asc = s1_asc.map(toNatural);
var s1_asc = s1_asc.map(RefinedLee);
var s1_asc = s1_asc.map(toDB);
var s1_asc_SD = s1_asc.reduce(ee.Reducer.stdDev());

var s1_dsc = ee.ImageCollection(s1_for_SD)
  .filterMetadata('orbitProperties_pass', 'equals' , 'DESCENDING');
var s1_dsc = s1_dsc.map(toGamma0);
var s1_dsc = s1_dsc.map(toNatural);
var s1_dsc = s1_dsc.map(RefinedLee);
var s1_dsc = s1_dsc.map(toDB);
var s1_dsc_SD = s1_dsc.reduce(ee.Reducer.stdDev());

var s1_SD = s1_asc_SD.add(s1_dsc_SD).divide(2);
var s1_SD = s1_SD.reproject(prj, null, 10);
var s1_SD_1 = s1_SD.clip(boreal_1);
var s1_SD_2 = s1_SD.clip(boreal_2);


var s1_mean =  ee.ImageCollection(s1_asc.merge(s1_dsc));
var s1_mean = s1_mean.mean();
var s1_mean = s1_mean.reproject(prj, null, 10);
var s1_mean_1 = s1_mean.clip(boreal_1);
var s1_mean_2 = s1_mean.clip(boreal_2);


var count_Boreal_50 = count_Boreal.reproject(prj, null, 50);
var count_Boreal_50 = count_Boreal_50.clip(boreal);

var count_Boreal_100 = count_Boreal.reproject(prj, null, 100);
var count_Boreal_100 = count_Boreal_100.clip(boreal);

///////////////////////////////////////////////////////////////////////////////////////
//2.6 Get standard deviation of pixel stack
//2.6 END
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
//2.7 Export images to Google Drive
//2.7 START
///////////////////////////////////////////////////////////////////////////////////////

Export.image.toDrive({
  image: s1_SD_1,
  description: 'S1_SD_boreal_1',
  scale: 10,
  region: boreal_1,
  maxPixels: 10E10
  
});
Export.image.toDrive({
  image: s1_SD_2,
  description: 'S1_SD_boreal_2',
  scale: 10,
  region: boreal_2,
  maxPixels: 10E10
  
});
Export.image.toDrive({
  image: PercentWater,
  description: 'HTV_Boreal',
  scale: 10,
  region: boreal,
  maxPixels: 10E10
  
});
Export.image.toDrive({
  image: s1_mean_1,
  description: 'S1_mean_boreal_1',
  scale: 10,
  region: boreal_1,
  maxPixels: 10E10
  
});

Export.image.toDrive({
  image: s1_mean_2,
  description: 'S1_mean_boreal_2',
  scale: 10,
  region: boreal_2,
  maxPixels: 10E10
  
});
Export.image.toDrive({
  image: count_Boreal_50,
  description: 'count_Boreal_50',
  scale: 50,
  region: boreal,
  maxPixels: 10E10
  
});
Export.image.toDrive({
  image: count_Boreal_100,
  description: 'count_Boreal_100',
  scale: 100,
  region: boreal,
  maxPixels: 10E10
  
});

///////////////////////////////////////////////////////////////////////////////////////
//2.7 Export images to Google Drive
//2.7 END
///////////////////////////////////////////////////////////////////////////////////////




//#####################################################################################
//#####################################################################################
///////////////////////////////////////////////////////////////////////////////////////
//2.Run Boreal HTV Algorithm
//2.END
///////////////////////////////////////////////////////////////////////////////////////
//#####################################################################################
//#####################################################################################









