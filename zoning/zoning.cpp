// zoning.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "classifier.h"

using namespace agrigate;

/*
input:
- csv
- .tiles
- params
- geojson

output:
- raster

implementation v1:
Algorithm
- select .tiles - v1 simple alg
- build RasterBuffer, nodata outside border
- classify
- postproc
- save raster

code:
- class agri::Classifier:
  - CreateNDVIBuffer (bbox,.tiles path)
  - ClipNDVIBuffer (geom, buffer)
  - Classify
  - SaveRaster


*/
//
//if -ranges == "", if o_gjson, than polygoni
//
//

/*
const list<GMXOptionDescriptor> listDescriptors = {
  { "-i", 0, 0, "input folder" },
  { "-o", 0, 0, "output folder" }

};


const list<string> listUsageExamples = {
  "ngw_lukoil -i c:\\CLS_input -o c:\\NGW_converted"
};
*/
const list<GMXOptionDescriptor> listDescriptors =
{
	//{ "-mean", 0, 0, "input mean csv" },
	//{ "-max", 0, 0, "input max csv" },
	{ "-v", 0, 0, "input geojson" },
	{"-sid",0, 0, "sceneid tiles file list"},
	{"-mosaic",1,0, "merge input ndvi tiles"},
	{ "-o_png", 0, 0, "output png" },
	{ "-o_tif", 0, 0, "output geotiff" },
	{"-m", 0, 0, "classify method: area, quantiles, intervals (deafault quantiles)"},
	{"-sig",0,0, "sigma parameter for equal intervals method (default 1.6)"},
	{"-o_geojson",0,0, "output vector zones"},
	{"-o_geojson_srs",0,0,"0-1 output vector zones srs code, (default 0 - epsg:4326)"},
	{ "-ncl", 0, 0, "number of classes (default 5)" },
	{"-filt", 0, 0, "majority filter iterations (default 2)"},
	{ "-z", 0, 0, "zoom (default = 13)" },
	{ "-ranges", 0, 0, "intervals to split" },
	{ "-col", 0, 1, "color table values" },
	{"-o_grid",0, 0, "output vector grid"},
	{ "-fids", 0, 0, "feature ids comma separated"}
};


const list<string> listUsageExamples =
{
	"zoning -v field.geojson -sid \"-S2A_20180811_080611_369_38UMC_ndvi.tiles,S2B_20180816_080559_274_38UMC_ndvi.tiles"
			" -o_png zones.png -o_png chveg_map.png -m area",
	"zoning -v field.geojson -sid \"NDVI/LC81720282014078LGN00_ndvi.tiles,NDVI/S2A_L1C_20160619_124_ndvi.tiles,"
	        "NDVI/LC81720282017102LGN00_ndvi.tiles\" -o_png zones.png -o_geojson zones.geojson -m area",
	"zoning -v field.geojson -sid \"NDVI/LC81720282014078LGN00_ndvi.tiles,NDVI/S2A_L1C_20160619_124_ndvi.tiles,"
			"NDVI/LC81720282017102LGN00_ndvi.tiles\" -o_png zones.png -o_geojson zones.geojson -m intervals -sig 2",

	//"zoning -mean mean.csv -max max.csv -v border.geojson -o_png zones.png -o_csv sceneid_list.csv",
	//"zoning -sid LC81720282014078LGN00,S2A_L1C_20160619_124,LC81720282017102LGN00 -v border.geojson -o_png zones.png -o_csv sceneid_list.csv",
	//"zoning -sid LC81720282014078LGN00,S2A_L1C_20160619_124,LC81720282017102LGN00 -v border.geojson -o_tif zones.tif -ranges 0.2,0.4,0.6,0.8,1"
	//zoning -v fields.shp -sid scene_id.csv -o_shp zones.shp -o_tif zones.tif -o_json zones.json
	//zoning -sid LC81720282014078LGN00,S2A_L1C_20160619_124,LC81720282017102LGN00 
	//-v border.geojson -o_tif zones.tif -ranges val1,val2,val3 -col 1 -col 2 -col 3 
	//-filt 2
	//classify method
	//zoning -ranges 101,111,121,131,136,141,146,151,156,161,166 -v "F:\mpotanin\data3\field_1.shp" -o_tif "F:\mpotanin\data3\field_1_zones.tif" -sid Y:\Kosmosnimki\Operative\Landsat8\ndvi\2017\2017-12-13\LC81750282017347LGN00_ndvi.tiles -z 12
};



#ifdef WIN32
//int main(int nArgs, wchar_t *pastrArgsW[])
int main(int nArgs, char* argv[])
{
	string *pastrArgs = new string[nArgs];

	for (int i = 0; i<nArgs; i++)
	{
		//GMXString::wstrToUtf8(pastrArgs[i], pastrArgsW[i]);
		pastrArgs[i] = argv[i];
		GMXString::ReplaceAll(pastrArgs[i], "\\", "/");
	}
#else
int main(int nArgs, char* argv[])
{
	string* pastrArgs = new string[nArgs];
	for (int i = 0; i<nArgs; i++)
		pastrArgs[i] = argv[i];
#endif

	if (!GMXGDALLoader::Load(GMXFileSys::GetPath(pastrArgs[0])))
	{
		cout << "ERROR: can't load GDAL" << endl;
		return 1;
	}

	cout << endl;

	if (nArgs == 1)
	{
		cout << "version: " << GMXFileSys::ReadTextFile(GMXFileSys::GetAbsolutePath(
			GMXFileSys::GetPath(pastrArgs[0]),
			"version.txt")
			) << endl;
		cout << "build date: " << __DATE__ << endl;
		GMXOptionParser::PrintUsage(listDescriptors, listUsageExamples);
		delete[]pastrArgs;
		return 0;
	}

/*
	Step 1: Parse cmd options, init main vars
*/

	GMXOptionParser oOptionParser;
	if (!oOptionParser.Init(listDescriptors, pastrArgs, nArgs))
	{
		cout << "ERROR: input cmd line is not valid" << endl;
		delete[]pastrArgs;
		return 2;
	}
	delete[]pastrArgs;


	double dblPixelOfset = 2;
	double dblSTDCoeff = oOptionParser.GetOptionValue("-sig") == "" ? 1.6 :
		atof(oOptionParser.GetOptionValue("-sig").c_str());
	int nWinSize = 5;
	int nClasses = oOptionParser.GetOptionValue("-ncl") == "" ? 5 : atoi(oOptionParser.GetOptionValue("-ncl").c_str());
	int nMajorityFilterNumRun = oOptionParser.GetOptionValue("-filt") == "" ? 0 :
								atoi(oOptionParser.GetOptionValue("-filt").c_str());

	list<string> listTileContainers;
	listTileContainers = GMXString::SplitCommaSeparatedText(oOptionParser.GetOptionValue("-sid"));
	if (listTileContainers.size() == 0)
	{
		cout << "ERROR: no image fits criterias" << endl;
		return 4;
	}

	int nZoom = oOptionParser.GetOptionValue("-z") == "" ? 0 : 
		atoi(oOptionParser.GetOptionValue("-z").c_str());


	string strClassifyMethod = oOptionParser.GetOptionValue("-m");



/*
	Step 2. Read input rasters (geotiff files or tile containers), init GeoRasterBuffer object
*/

	OGRGeometry* poVecBorder = 0;
	GeoRasterBuffer* poGeoBuffer = 0;
	int nFIDCount = 0;
	int* panFIDs = 0;
	if (oOptionParser.GetOptionValue("-fids") != "")
	{
		list<string> listFIDs = GMXString::SplitCommaSeparatedText(oOptionParser.GetOptionValue("-fids"));
		panFIDs = new int[listFIDs.size()];
		for (auto strFID : listFIDs)
		{
			panFIDs[nFIDCount] = atoi(strFID.c_str());
			nFIDCount++;
		}
	}

	if (((*listTileContainers.begin()).find(".tiles") != string::npos) ||
		((*listTileContainers.begin()).find(".gmxtiles") != string::npos))
	{

		gmx::MercatorTileMatrixSet oMercTMS;
		poVecBorder = gmx::VectorOperations::ReadIntoSingleMultiPolygon(
			oOptionParser.GetOptionValue("-v"), oMercTMS.GetTilingSRSRef(), panFIDs, nFIDCount);

		poGeoBuffer = GeoRasterBuffer::InitFromNDVITiles(listTileContainers,
														poVecBorder,
														strClassifyMethod == "quantiles" ? false : true,
														oOptionParser.GetOptionValue("-mosaic") != "",
														-0.5, nZoom);

		if (!poGeoBuffer)
		{
			cout << "ERROR: InitFromNDVITiles\n";
			return 5;
		}
	}
	else
	{
		poGeoBuffer = GeoRasterBuffer::InitFromNDVIFiles(listTileContainers,
			oOptionParser.GetOptionValue("-v"),
			false,
			-0.5);

		if (!poGeoBuffer)
		{
			cout << "ERROR: InitFromNDVIFiles\n";
			return 5;
		}

		poVecBorder = gmx::VectorOperations::ReadIntoSingleMultiPolygon(
						oOptionParser.GetOptionValue("-v"), poGeoBuffer->GetSRSRef(), panFIDs, nFIDCount);
	}
	
	
	
/*
	Step 3. Classify raster data, create ClassifiedRasterBuffer object
*/

	int* panIntervals = 0;
	ClassifiedRasterBuffer* poClassifiedBuffer = 0;
	if (strClassifyMethod == "quantiles")
	{
		nClasses = 5;
		double padblQuantiles[6] = {0,0.25,0.4,0.6,0.75,1.};
		poClassifiedBuffer = poGeoBuffer->ClassifyByPredefinedQuantiles(nClasses, 
																		padblQuantiles);
	}
	else if (strClassifyMethod == "intervals")
		poClassifiedBuffer = poGeoBuffer->ClassifyEqualIntervals(nClasses, panIntervals, dblSTDCoeff);
	else if (strClassifyMethod == "area" || strClassifyMethod == "")
		poClassifiedBuffer = poGeoBuffer->ClassifyEqualArea(nClasses, panIntervals);
	else if (strClassifyMethod == "raw")
	{
		panIntervals = new int[101];
		nClasses = 100;
		for (int i=0;i<=100;i++)
			panIntervals[i]=i;
		poClassifiedBuffer = poGeoBuffer->ClassifyByPredefinedIntervals(nClasses, panIntervals);
	}
	else
	{
		list<string> listInterval = GMXString::SplitCommaSeparatedText(strClassifyMethod);
		
		panIntervals = new int[listInterval.size()+2];
		nClasses = listInterval.size() + 1;
		int i = 1;
		panIntervals[0] = 0;
		panIntervals[nClasses] = 100*poGeoBuffer->get_num_bands();
		for (auto str : listInterval)
		{
			panIntervals[i] = (int)(100 * atof(str.c_str()) +0.5);
			if (panIntervals[i] <= 0 && panIntervals[i] > 99)
			{
				cout << "ERROR: -m parameter is not valid: " << strClassifyMethod << endl;
				return 2;
			}
			i++;
		}

		if (panIntervals[2] <= panIntervals[1])
		{
			cout << "ERROR: -m parameter is not valid: " << strClassifyMethod << endl;
			return 2;
		}
		poClassifiedBuffer = poGeoBuffer->ClassifyByPredefinedIntervals(nClasses, panIntervals);
	}
	
	if (poClassifiedBuffer == 0)
	{
		cout << "ERROR: classify error" << endl;
		return 5;
	}
  

/*
	Step 4. Postprocess classification results: 
	   - interpolate pixels along vector border
	   - apply majority filter
	   - adjust extent
*/

	poClassifiedBuffer->ReplaceByInterpolatedValues(poVecBorder, 1, 1);
	
	for (int i = 0; i < nMajorityFilterNumRun; i++)
		poClassifiedBuffer->ApplyMajorityFilter(nWinSize);
	  
	poClassifiedBuffer->AdjustExtentToClippedArea();
	
  
/*
	Step 4. Postprocess classification results:
		- interpolate pixels along vector border
		- apply majority filter
		- adjust extent
*/

	if (oOptionParser.GetOptionValue("-o_geojson") != "")
	{

		map<int, OGRMultiPolygon*> mapZones = poClassifiedBuffer->Polygonize();
		
		ZoningMap oZM;
		oZM.InitDirectly(mapZones);
		//debug
		//oZM.SaveToVectorFile(oOptionParser.GetOptionValue("-o_geojson"), poClassifiedBuffer->GetSRSRef());
		//end-debug


		oZM.Clip(poVecBorder);  //should fix it
		
		//OGRSpatialReference oSRS;
		OGRSpatialReference oWGS84SRS;
		
		oWGS84SRS.SetWellKnownGeogCS("EPSG:4326");
		//oSRS.importFromProj4("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs");
		//oZM.TransformToSRS(oMercTMS.GetTilingSRSRef(), &oSRS);
		//char* pachWKT;
		//poClassifiedBuffer->GetSRSRef()->exportToWkt(&pachWKT);

		//oZM.TransformToSRS(poClassifiedBuffer->GetSRSRef(), &oSRS);
		//oZM.SaveToVectorFile(oOptionParser.GetOptionValue("-o_geojson"),&oSRS);
		oZM.SaveToVectorFile(oOptionParser.GetOptionValue("-o_geojson"), poClassifiedBuffer->GetSRSRef());
		//oZM.TransformToSRS(poClassifiedBuffer->GetSRSRef(), &oWGS84SRS);
		//oZM.SaveToVectorFile(oOptionParser.GetOptionValue("-o_geojson"),&oWGS84SRS);

		//OGRSpatialReference oSRS;
		//oWGS84SRS.SetWellKnownGeogCS("EPSG:4326");
		//oZM.TransformToSRS(oMercTMS.GetTilingSRSRef(), &oWGS84SRS);
		//oZM.SaveToVectorFile(oOptionParser.GetOptionValue("-o_geojson"),&oWGS84SRS);

		//oZM.SaveToFile(oOptionParser.GetOptionValue("-o_geojson"));
		//oZM.FilterByArea(5000);
		//oZM.SaveToFile("F:\\mpotanin\\byukreevka\\zones3_clipped_nofilt_field6.shp", "");
	}


/*
	Step 5. Save zones in raster fromat, apply color palette  
*/
	string strRasterOutput = "";
	if (oOptionParser.GetOptionValue("-o_png") != "")
		strRasterOutput = oOptionParser.GetOptionValue("-o_png");
	else if (oOptionParser.GetOptionValue("-o_tif")!="")
		strRasterOutput = oOptionParser.GetOptionValue("-o_tif");
	
	if (strRasterOutput != "")
	{
		GDALColorTable *poColTab = new GDALColorTable(GPI_RGB);

		if ((oOptionParser.GetValueList("-col").size() == 0) && (oOptionParser.GetOptionValue("-ranges") == ""))
		{
			if (nClasses == 3)
			{
				GDALColorEntry arrColors[6];
				arrColors[0].c1 = 0; arrColors[0].c2 = 0; arrColors[0].c3 = 0;
				arrColors[1].c1 = 255; arrColors[1].c2 = 0; arrColors[1].c3 = 0;
				arrColors[2].c1 = 212; arrColors[2].c2 = 255; arrColors[2].c3 = 190;
				arrColors[3].c1 = 47; arrColors[3].c2 = 140; arrColors[3].c3 = 30;
				for (int i = 0; i <= nClasses; i++)
					poColTab->SetColorEntry(i, &arrColors[i]);
			}
			else
			{
				GDALColorEntry arrColors[6];
				arrColors[0].c1 = 0; arrColors[0].c2 = 0; arrColors[0].c3 = 0;
				arrColors[1].c1 = 255; arrColors[1].c2 = 0; arrColors[1].c3 = 0;
				arrColors[2].c1 = 247; arrColors[2].c2 = 209; arrColors[2].c3 = 59;
				arrColors[3].c1 = 212; arrColors[3].c2 = 255; arrColors[3].c3 = 190;
				arrColors[4].c1 = 76; arrColors[4].c2 = 227; arrColors[4].c3 = 0;
				arrColors[5].c1 = 47; arrColors[5].c2 = 140; arrColors[5].c3 = 30;
				for (int i = 0; i <= nClasses; i++)
					poColTab->SetColorEntry(i, &arrColors[i]);
			}

		}
		else
		{
			GDALColorEntry oColor;
			oColor.c1 = 0; oColor.c2 = 0; oColor.c3 = 0;
			poColTab->SetColorEntry(0, &oColor);

			int i = 0;
			unsigned char panRGB[3];
			for (string strColor : oOptionParser.GetValueList("-col"))
			{
				GMXString::ConvertStringToRGB(strColor, panRGB);
				oColor.c1 = panRGB[0]; oColor.c2 = panRGB[1]; oColor.c3 = panRGB[2];
				poColTab->SetColorEntry(i + 1, &oColor);
				i++;
			}
		}
		
		poClassifiedBuffer->set_color_table(poColTab);

		poClassifiedBuffer->SaveGeoRefFile(strRasterOutput);
		delete(poColTab);
	}

/*
	Step 6 (optional). Save grid file
*/

	if (oOptionParser.GetOptionValue("-o_grid") != "")
	{
		OGRFeature** paoFeatures = 0;
		int nFeatures = VectorOperations::ReadAllFeatures(oOptionParser.GetOptionValue("-v"), 
															paoFeatures, 
															poClassifiedBuffer->GetSRSRef());

		double* padblClassCount = new double[poClassifiedBuffer->GetNumClasses()];
		for (int i = 0; i<poClassifiedBuffer->GetNumClasses(); i++)
			padblClassCount[i]=0.;

		for (int i = 0; i < nFeatures; i++)
		{
			if (panFIDs)
			{
			int j = 0;
			for (j=0;j<nFIDCount;j++)
				if (panFIDs[j]==i+1) break;
			if (j==nFIDCount) continue;
			}

			string strGridFile = GMXFileSys::GetAbsolutePath(oOptionParser.GetOptionValue("-o_grid"),
															"grid_field_" + GMXString::ConvertIntToString(i+1) + ".shp");
			poClassifiedBuffer->PolygonizePixels(strGridFile, paoFeatures[i]->GetGeometryRef(), true);
			OGRFeature::DestroyFeature(paoFeatures[i]);

			OGRFeature** paoCellFeatures;
			int nCells = VectorOperations::ReadAllFeatures(strGridFile,paoCellFeatures);
			for (int j = 0; j < nCells; j++)
			{
			padblClassCount[paoCellFeatures[j]->GetFieldAsInteger("DN")-1]+=paoCellFeatures[j]->GetFieldAsDouble("FRACTION");
			delete(paoCellFeatures[j]);
			}
			delete[]paoCellFeatures;
		}
		string strClassesCountFile = GMXFileSys::GetAbsolutePath(oOptionParser.GetOptionValue("-o_grid"),"classes_count.txt");
		FILE* fp = fopen(strClassesCountFile.c_str(), "w");
		for (int i = 0; i<poClassifiedBuffer->GetNumClasses();i++)
		  fprintf(fp, "%.2lf\n", padblClassCount[i]);
		fclose(fp);

		delete[]paoFeatures;
		delete[]padblClassCount;
	}


	delete[]panFIDs;
	delete[]panIntervals;
	delete(poClassifiedBuffer);
	delete(poGeoBuffer);



	return 0;
}

