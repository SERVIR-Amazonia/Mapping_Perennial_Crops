//Importar os arquivos
var Area_Geral = ee.FeatureCollection("users/tassioigawa/Curso_Culturas_Perenes/AI_mun");
var Pol_Dende = ee.FeatureCollection("users/tassioigawa/Curso_Culturas_Perenes/Pol_Dende");
var Pol_N_Dende = ee.FeatureCollection("users/tassioigawa/Curso_Culturas_Perenes/Pol_N_Dende");
var Pts_Dende = ee.FeatureCollection("users/tassioigawa/Curso_Culturas_Perenes/Pts_Dende");
var Pts_N_Dende = ee.FeatureCollection("users/tassioigawa/Curso_Culturas_Perenes/Pts_N_Dende");

//Área de estudo
var area_estudo = Area_Geral;

//Amostras
var pts_dende = Pts_Dende;
var pts_ndende= Pts_N_Dende;
var dende = Pol_Dende;
var ndende = Pol_N_Dende;
var features = Pol_Dende.merge(Pol_N_Dende);
var Amostrasdende = Pts_Dende.merge(Pts_N_Dende);
print('Amostras', features)

//Carregar ALOS
var alos = ee.ImageCollection("JAXA/ALOS/PALSAR/YEARLY/SAR");
//Carregar Sentinel-1 TOA reflectance data/
var s1 = ee.ImageCollection("COPERNICUS/S1_GRD");
//Carregar Sentinel-2 TOA reflectance data/
var s2 = ee.ImageCollection('COPERNICUS/S2');
//Carregar Landsat 8 TOA reflectance data/
var lsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA');

///////////////////////////////////////////////// Funções inserção de variaveis ////////////////////////////////////
//Função para rodar o filtro lee (ALOS 2 -PALSAR-2)
function filterleealos(image)
{
  // Filtro Lee para suavizar e reduzir os ruídos das imagens de Radar (Homogenius Class - Enhanced Lee Filter)
  var kernellee = ee.Kernel.square(3,'pixels');
  var BandsLeealos =['HHLee','HVLee']
  var leefiltersalos = image.reduceNeighborhood({
    reducer:ee.Reducer.mean(),
    kernel:kernellee,
  });
  return image.addBands(leefiltersalos.rename(BandsLeealos));
}

//Função para rodar o filtro lee (Sentinel 1)
function filterlee(image)
{
  // Filtro Lee para suavizar e reduzir os ruídos das imagens de Radar (Classe Homogênea - Filtro Lee Melhorado)
  var kernellee = ee.Kernel.square(3,'pixels');
  var BandsLee =['VVLee','VHLee']
  var leefilters = image.reduceNeighborhood({
    reducer:ee.Reducer.mean(),
    kernel:kernellee,
  });
  return image.addBands(leefilters.rename(BandsLee));
}


//Função para Adicionar o RVI (ALOS 2-PALSAR-2)
function addalosRVI(image)
{
  var rvialos = image.expression(
  '4 * hv / (hv + hh)',
  {
     hh: image.select('HHLee'),    
     hv: image.select('HVLee')
  });
  return image.addBands(rvialos.rename('RVIAlos'));
}
 
//Função para Adicionar o RVI (Sentinel 1)
function addS1RVI(image)
{
  var rvi = image.expression(
  '4 * vh / (vh + vv)',
  {
     vv: image.select('VVLee'),    
     vh: image.select('VHLee')
  });
  return image.addBands(rvi.rename('RVI'));
}


//Função para Adicionar o NDVI do Landsat 8
function addLandsatNDVI(image)
{
  var Landsat_NDVI = image.expression(
    '((ir-red)/(ir+red))',
    {
      ir: image.select('B5'),
      red: image.select('B4')
    });
    return image.addBands(Landsat_NDVI.rename('LS8_NDVI'));
}


//LAI
function addLAI(image)
{
    var ndvi = image.expression(
    '((red - ir) / (red + ir))',
    {
        red: image.select('B8'),    
        ir: image.select('B4')});
        
    var LAI = ndvi.expression('0.57 * exp(2.33 * b(0))').rename('LAI');

    image = image.addBands(LAI.rename('LAI'));
return image;
}


// Função para a criação das máscara de nuvem usando Sentinel-2 QA60 band.
function maskS2clouds(image) {
  var qa = image.select('QA60');

// Bits 10 and 11 are nuvens and cirrus, respectivamente.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  
// Ambos os sinalizadores devem ser definidos como zero, indicando condições claras.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}


//função para adicionar indices de vegetação
function addIndices(image){
  return image.addBands([
    image.select('B8','B4').normalizedDifference().multiply(10000).add(10000).rename('NDVI')]);
}

                     
// Fazer o empilhamento de imagens
var alospalsar =  alos.filterDate('2017-06-01', '2018-10-31')
                  .select(['HH', 'HV'])
                  .map(filterleealos)
                  .map(addalosRVI);

var landsat8 =  lsat8.filterDate('2020-06-01', '2020-10-31')
                  .filterMetadata('CLOUD_COVER','less_than', 25)  
                  .select(['B4', 'B5'])
                  .map(addLandsatNDVI);

var sentinel1 =  s1.filterDate('2020-06-01', '2020-10-31')
                  .select(['VV', 'VH'])
                  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                  .filter(ee.Filter.eq('instrumentMode', 'IW'))
                  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                  .map(filterlee)
                  .map(addS1RVI);
                  
var sentinel2 = s2.filterDate('2020-06-01', '2020-10-31')
                  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than', 10)  
                  .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6','B7','B8', 'B8A', 'B11', 'B12','QA60'])
                  .map(maskS2clouds)
                  .map(addIndices)
                  .map(addLAI);

// Composição de um conjunto de dados na área de estudo

var landsat8 = landsat8.median().clip(area_estudo);
var s1lee = sentinel1.median().clip(area_estudo);
var s2filtered = sentinel2.median().clip(area_estudo);
var aloslee = alospalsar.median().clip(area_estudo);
var stack = s2filtered.addBands(landsat8).addBands(s1lee).addBands(aloslee);

// Mostrar o stack de imagens
Map.addLayer(stack, {bands: ['B11', 'B8', 'B4'], min: 0, max: 0.3}, 'stack', false);
//Map.addLayer(al2_composite, {bands: ['HH', 'HV'], min: 0, max: 8500}, 'ALOS2 composite', false);

var bands = ['B1', 'B2', 'B3','B4', 'B5', 'B6','B7','B8','B8A', 'B11', 'B12','LAI','LS8_NDVI','VV','VH','RVI','HH','HV','RVIAlos'];
var classe = 'Classe';
var presence = 'presence';

////////////////////////////////////////// Amostras /////////////////////////////////////////////////////
var training = stack.select(bands).sampleRegions({
  collection: features,
  properties: [classe],
  scale: 50
});

var trainingmaxent = stack.sampleRegions({
  collection: Amostrasdende,
  properties: [presence],
  scale: 50
});

//Acurária 
// Opcionalmente, faça alguma avaliação de precisão. Adicione uma coluna de
// uniformes aleatórios para o conjunto de dados de treinamento.
var withRandom = training.randomColumn('random');

//Queremos reservar alguns dos dados para teste, para evitar overfitting do modelo.
var split = 0.7;  // Aproximadamente 70% treinando, 30% testando.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));


//////////////////////////////////////// Treinamentos dos Classificadores //////////////////////////////////////////

// Treinado com 70% das amostras
var trainedClassifier = ee.Classifier.smileRandomForest(100).train({
  features: trainingPartition,
  classProperty: classe,
  inputProperties: bands
});


// Definir e treinar o classificador Maxent de uma imagem por amostragens de pontos.
var classifier = ee.Classifier.amnhMaxent().train({
  features: trainingmaxent,
  classProperty: 'presence',
  inputProperties: stack.bandNames()
});

//Treinando o Classificador (Classification and Regression Trees - CART)
var trained = ee.Classifier.smileCart().train(trainingPartition, classe, bands);

////////////////////////////////////////////// Imprimir os resultados /////////////////////////////////////////////

//Imprima algumas informações sobre o RF
print('Random Forest, explicado', trainedClassifier.explain());

//Imprima algumas informações sobre o CART
print('CART, explicado', trained.explain());

//Imprima algumas informações sobre o MaxEnt
print('MaxEnt, explicado', classifier.explain());

//////////////////////////////////////////// Classificadores //////////////////////////////////////////////////////

//Classifique o composto (Random Forest)
var classificado = stack.classify(trainedClassifier);

//Classifique o composto (CART)
var classificado_cart = stack.select(bands).classify(trained);

// Classificar o composto (MaxEnt)
var imageClassifiedMax = stack.classify(classifier);

// Classifica o teste na Coleção de Feições (Random Forest)
var test = testingPartition.classify(trainedClassifier);

// Classifica o teste na Coleção de Feições (CART)
var testCart = testingPartition.classify(trained);

////////////////////////////////////////// Visualizar Camadas ///////////////////////////////////////////////////

// Adicionando Camadas 
Map.centerObject(area_estudo,8)
Map.addLayer(classificado_cart, {min: 1, max: 2, palette: ['green','blue']},'CART');
Map.addLayer(classificado, {min: 1, max: 2, palette: ['green','blue']},'Random Forest');
Map.addLayer(imageClassifiedMax, {bands: 'probability', min: 0, max: 1}, 'Probabilidade - MaxEnt');


// Definir zonas do modelo probabilístico MaxENT
var MaxEnt = imageClassifiedMax.select('probability');
var zones = MaxEnt.gt(0.01).add(MaxEnt.gt(0.616));
zones = zones.updateMask(zones.neq(0));


//Visualizar as Zonas do MaxENT
Map.addLayer(zones, {min: 1, max: 2, palette: ['0000FF', '00FF00']}, 'Zonas - MaxEnt');

//////////////////////////////////////////Filtro Majoritário na Classificação////////////////////////////

// CART
var kernel = ee.Kernel.square(2);
var fmajoritclasscart = classificado_cart.reduceNeighborhood({
  reducer:ee.Reducer.mode(),
  kernel:kernel,
});

// MaxEnt
var fmajoritclassmax = zones.reduceNeighborhood({
  reducer:ee.Reducer.mode(),
  kernel:kernel,
});

// Random Forest
var fmajoritclassrf = classificado.reduceNeighborhood({
  reducer:ee.Reducer.mode(),
  kernel:kernel,
});

//Adicionando a camada
Map.addLayer(fmajoritclasscart,{min: 1, max: 2, palette: ['black', 'green']},'Majoritario_CART');
Map.addLayer(fmajoritclassmax,{min: 1, max: 2, palette: ['black', 'green']},'Majoritario_Max');
Map.addLayer(fmajoritclassrf,{min: 1, max: 2, palette: ['black', 'green']},'Majoritario_RF');


//////////////////////////////////////////// Validação /////////////////////////////////////////////////////
// MaxEnt
var verdpos = zones.eq(2).clip(dende);
var falpos = zones.eq(2).clip(ndende);
var verdneg = zones.eq(1).clip(dende);
var falng = zones.eq(1).clip(ndende);
var areatotaldende = zones.gte(1).clip(dende);
var areatotalndende = zones.gte(1).clip(ndende);

var verdposarea = verdpos.multiply(ee.Image.pixelArea());
var verdnegarea = verdneg.multiply(ee.Image.pixelArea());
var falposarea = falpos.multiply(ee.Image.pixelArea());
var falngarea = falng.multiply(ee.Image.pixelArea());
var areatotalareadende = areatotaldende.multiply(ee.Image.pixelArea());
var areatotalareandende = areatotalndende.multiply(ee.Image.pixelArea());


var VParea = verdposarea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: dende.geometry(),
  scale: 10,
  maxPixels: 1e13
  });

var VNarea = verdnegarea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: dende.geometry(),
  scale: 10,
  maxPixels: 1e13
  });
  
var FNarea = falngarea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ndende.geometry(),
  scale: 10,
  maxPixels: 1e13
  });

var FParea = falposarea.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ndende.geometry(),
  scale: 10,
  maxPixels: 1e13
  });
  
var ATareadende = areatotalareadende.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: dende.geometry(),
  scale: 10,
  maxPixels: 1e13
  }); 

var ATareandende = areatotalareandende.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ndende.geometry(),
  scale: 10,
  maxPixels: 1e13
  }); 


var VPSqKm = ee.Number(
  VParea.get('probability'));

var VNSqKm = ee.Number(
  VNarea.get('probability'));
  
var FPSqKm = ee.Number(
  FNarea.get('probability'));
  
var FNSqKm = ee.Number(
  FParea.get('probability'));
  
var ATdendeSqKm = ee.Number(
  ATareadende.get('probability'));

var ATndendeSqKm = ee.Number(
  ATareandende.get('probability'));
  
var acuracia = ee.Number(VPSqKm).divide(ee.Number(ATdendeSqKm))  
var areatotal = ee.Number(ATdendeSqKm).add(ee.Number(ATndendeSqKm))

//Tx de Aceitação Relativa
var po = (ee.Number(VPSqKm).add(ee.Number(FNSqKm))).divide(ee.Number(areatotal));
//Probabilidade de ambos randomicamente aceitarem um projeto
var pac = ((ee.Number(VPSqKm).add(ee.Number(FPSqKm))).divide(ee.Number(areatotal))).multiply((ee.Number(VPSqKm).add(ee.Number(VNSqKm))).divide(ee.Number(areatotal)));
//Probabilidade de ambos randomicamente rejeitarem
var prej = ((ee.Number(FNSqKm).add(ee.Number(VNSqKm))).divide(ee.Number(areatotal))).multiply((ee.Number(FPSqKm).add(ee.Number(FNSqKm))).divide(ee.Number(areatotal)));
//Taxa hipotética de aceitação
var pe = ee.Number(pac).add(ee.Number(prej));
//kappa - Nesse caso, não apresentou concordância devido o valor abaixo de 0
var kappaMax = (ee.Number(po).subtract(ee.Number(pe))).divide(ee.Number(1).subtract(ee.Number(pe)));

print('Acurácia geral da validação - MaxEnt', acuracia);
print('Kappa - MaxEnt', kappaMax);

//Random Forest
var MatrizErro = test.errorMatrix(classe, 'classification');
var AG = MatrizErro.accuracy()
var Kappa = MatrizErro.kappa()

print('Acurácia geral da validação - RF', AG)
print('Kappa - RF',Kappa)

//Classification and Regression Trees
var MatrizErroCart = testCart.errorMatrix(classe, 'classification');
var AG = MatrizErroCart.accuracy()
var Kappa = MatrizErroCart.kappa()

print('Acurácia geral da validação - CART', AG)
print('Kappa - CART',Kappa)

///////////////////////////////////////Visualização dos Modelos////////////////////////////////////////////
//Código - Fonte: Earth Engine Apps
// Cria uma variável com as coordenadas do ponto
var pt1 = /* color: #ffff00 */ee.Geometry.Point([-48.1933, -1.3658]);

var coordenadas = pt1.geometries();
var Vale = coordenadas.get(0);
print('Córrego do Feijão', Vale);

//Selecionar as imagens 
var images=[
  ee.ImageCollection(classificado),
  ee.ImageCollection(classificado_cart),
  ee.ImageCollection(zones),
  ee.ImageCollection(stack)
  ];

var vis = [
  {min: 1, max: 2, palette: ['#000000','#DC143C']},
  {min: 1, max: 2, palette: ['#000000','#DC143C']},
  {min: 1, max: 2, palette: ['#000000','#DC143C']},
  {min: 0, max: 0.3, bands: ['B4', 'B3', 'B2']}
];

var NAMES = [
  'Random Forest',
  'CART',
  'MaxEnt',
  'Imagens'
];

// Cria um mapa para cada opção de visualização.
var maps = [];
NAMES.forEach(function(name,index) {
  var map = ui.Map();
  map.add(ui.Label(name));
  map.addLayer(images[index], vis[index], name);
  map.setControlVisibility(true);
  maps.push(map);
});

var linker = ui.Map.Linker(maps);
// Habilita o Zoom no mapa da esquerda "Antes".
maps[0].setControlVisibility({zoomControl: true});
// Mostra a escala (por exemplo, 500m) no mapa da direita.
maps[0].setControlVisibility({scaleControl: true});
maps[1].setControlVisibility({scaleControl: true});
maps[2].setControlVisibility({scaleControl: true});
maps[3].setControlVisibility({scaleControl: true});

// Cria o título.
var title = ui.Label('Modelos - Dendê', {
  stretch: 'horizontal',
  textAlign: 'center',
  fontWeight: 'bold',
  fontSize: '24px',
  color: '#0000FF'
});

// Cria um painel para armazenar a interface gráfica.
var panel = ui.Panel();
panel.style().set('width', '300px');

// Cria uma grade de mapas.
var mapGrid = ui.Panel([
    ui.Panel([maps[0], maps[1]], null, {stretch: 'both'}),
    ui.Panel([maps[2], maps[3]], null, {stretch: 'both'})
  ],
  ui.Panel.Layout.Flow('horizontal'), {stretch: 'both'}
);

// Adiciona os mapas e títulos na "ui.root".
ui.root.widgets().reset([title, mapGrid]);
ui.root.setLayout(ui.Panel.Layout.Flow('vertical'));

// Centraliza os mapas.
maps[0].setCenter(-48.8828, -2.553,12);

ui.Map.Linker([maps[0], maps[1], maps[2]]);

//////////////////////////////////////////// Exportação das Classificações //////////////////////////////////
//CART
Export.image.toDrive({
  image:fmajoritclasscart,
  description:'Dende_CART',
  scale:10,
  folder:'Download_gee',
  region:area_estudo,
  maxPixels:1e13
});

// MaxEnt
Export.image.toDrive({
  image:fmajoritclassmax,
  description:'Dende_MaxEnt',
  scale:10,
  folder:'Download_gee',
  region:area_estudo,
  maxPixels:1e13
});

//Random forest 
Export.image.toDrive({
  image:fmajoritclassrf,
  description:'Dende_RF',
  scale:10,
  folder:'GEE',
  region:area_estudo,
  maxPixels:1e13
});

//Amostras
Export.table.toDrive({
  collection: dende,
  description:'Dende_Amostras_pol',
  folder: 'GEE',
  fileFormat: 'KML'
});

