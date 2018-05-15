/////////////////////////////////////////////////////////////////////////////////////////////
//                         HarvestAreaRecovery                                             //
//                                                                                         //
//Description: Harvest area recovery for all harvest areas in the ABMI HFI 2016            //
//             measured with Landsat 5 time series                                         //
//                                                                                         //
//Author: Evan R. DeLancey, Alberta biodiversity monitoring insitute                       //
//        Jen Hird, University of Calgary                                                  //
/////////////////////////////////////////////////////////////////////////////////////////////

//filter harvest areas with orgin year over 1985
var HA = table.filterMetadata('YEAR', 'less_than', 1999);



//add harvest area polygons to map
Map.addLayer(HA, {}, 'Harvest Areas');



//get Landsat 5 for All of Alberta
var L5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
  .filterBounds(geometry);
  


//List of year to plot
var years = ee.List.sequence(1985, 2011);

//function for clipping images by chosen harvest area
var Clip = function(image) {
  return image.clip(geometry);
};

var L5 = L5.map(Clip);

//take only "clear" pixels
var maskcloud = function(image) {
  var diff = image.select(['sr_cloud_qa']);
  return image.updateMask(diff.lt(2));
};
var L5 = L5.map(maskcloud);


// Function to calculate Tasseled Cap Brightness (based on Crist, 1985)
var addTCB = function(image) {
  var TCB = image.expression(
    '(0.2043*B1) + (0.4158*B2) + (0.5524*B3) + (0.5741*B4) + (0.3124*B5) + (0.2303*B7)',{
      'B1' : image.select(['B1']),
      'B2' : image.select(['B2']),
      'B3' : image.select(['B3']),
      'B4' : image.select(['B4']),
      'B5' : image.select(['B5']),
      'B7' : image.select(['B7'])
    }).rename('TCB');
  return image.addBands(TCB);
};
Series = L5.map(addTCB);


//function for computing median TCB per summer season
var seasonalMean = years.map(function(y) {
  var start = ee.Date.fromYMD(y, 5, 15);
  var stop = start.advance(4, 'month');
  var median = Series.filterDate(start, stop).median();
  return median.set('year', y);
});

var Series = ee.ImageCollection(seasonalMean);
print(Series);
var Albedo_Series = Series.select(['TCB']);

var subset = HA.filterBounds(geometry);

//Build table with median TCB for every year, for every harvest area
var triplets_Albedo = Albedo_Series.map(function(image) {
  return image.select('TCB').reduceRegions({
    collection: subset.select(['OBJECTID']), 
    reducer: ee.Reducer.mean(), 
    scale: 30
  }).filter(ee.Filter.neq('mean', null))
    .map(function(f) { 
      return f.set('year', image.get('year'));
    });
}).flatten();


//Export csv to drive to plot in ggplot2
Export.table.toDrive({
  collection: triplets_Albedo, 
  description: 'HAR_DF_1986_1998', 
  fileFormat: 'CSV'
});


//DONE SCRIPT
