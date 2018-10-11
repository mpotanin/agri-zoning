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
int nDescriptors = 14;
const GMXOptionDescriptor asDescriptors[] =
{
	//{ "-mean", 0, 0, "input mean csv" },
	//{ "-max", 0, 0, "input max csv" },
	{ "-v", 0, 0, "input geojson" },
	{"-sid",0, 0, "sceneid tiles file list"},
	{ "-o_png", 0, 0, "output png" },
	{ "-o_tif", 0, 0, "output geotiff" },
	//{ "-o_csv", 0, 0, "output csv file" },
	{"-m", 0, 0, "classify method: area, intervals (deafault area)"},
	{"-sig",0,0, "sigma parameter for equal intervals method (default 1.6)"},
	{"-o_geojson",0,0, "output vector zones"},
	{"-o_geojson_srs",0,0,"0-1 output vector zones srs code, (default 0 - epsg:4326)"},
	{ "-ncl", 0, 0, "number of classes (default 5)" },
	{"-filt", 0, 0, "majority filter iterations (default 2)"},
	{ "-z", 0, 0, "zoom (default = 13)" },
	{ "-ranges", 0, 0, "intervals to split" },
	{ "-col", 0, 1, "color table values" },
	{"-o_p2p",0, 0, "output vector file for pixels to polygons operation"}
};

int nExamples = 2;
const string astrUsageExamples[] =
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
int _tmain(int nArgs, wchar_t *pastrArgsW[])
{


	string *pastrArgs = new string[nArgs];
	for (int i = 0; i<nArgs; i++)
	{
		GMXString::wstrToUtf8(pastrArgs[i], pastrArgsW[i]);
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
		GMXOptionParser::PrintUsage(asDescriptors, nDescriptors, astrUsageExamples, nExamples);
		delete[]pastrArgs;
		return 0;
	}

	GMXOptionParser oOptionParser;
	if (!oOptionParser.Init(asDescriptors, nDescriptors, pastrArgs, nArgs))
	{
		cout << "ERROR: input cmd line is not valid" << endl;
		delete[]pastrArgs;
		return 2;
	}
	delete[]pastrArgs;


	////////////////////////////////////////////////////////////////////////////
	string strAliasPath = "kosmosnimki";
	string strRealPath = "//tinkerbell-smb/ifs/kosmosnimki";

	double dblPixelOfset = 2;
	double dblSTDCoeff = oOptionParser.GetOptionValue("-sig") == "" ? 1.6 :
		atof(oOptionParser.GetOptionValue("-sig").c_str());
	int nWinSize = 5;
	int nClasses = oOptionParser.GetOptionValue("-ncl") == "" ? 5 : atoi(oOptionParser.GetOptionValue("-ncl").c_str());
	int nMajorityFilterNumRun = oOptionParser.GetOptionValue("-filt") == "" ? 2 :
								atoi(oOptionParser.GetOptionValue("-filt").c_str());
	
	////////////////////////////////////////////////////////////////////////////


	list<string> listTileContainers;
	if (oOptionParser.GetOptionValue("-sid") == "")
	{
		agrigate::NDVIProfile oNDVIProfile;
		if (!oNDVIProfile.ParseInputData(oOptionParser.GetOptionValue("-mean"),
			oOptionParser.GetOptionValue("-max"),
			strAliasPath,
			strRealPath))
		{
			cout << "ERROR: can't parse input mean, max csv files" << endl;
			return 3;
		}

		listTileContainers = oNDVIProfile.SelectInputForClassification();
	}
	else
	{
		listTileContainers = GMXString::SplitCommaSeparatedText(oOptionParser.GetOptionValue("-sid"));
		if (strAliasPath != "")
		{
			for (string iter : listTileContainers)
			{
				if (iter.find(strAliasPath) == 0)
					iter = strRealPath + iter.substr(strAliasPath.size());
			}
		}
	}
		
	if (listTileContainers.size() == 0)
	{
		cout << "ERROR: no image fits criterias" << endl;
		return 4;
	}

	int nZoom = oOptionParser.GetOptionValue("-z") == "" ? 0 : 
		atoi(oOptionParser.GetOptionValue("-z").c_str());

	gmx::MercatorTileMatrixSet oMercTMS;


	OGRGeometry* poVecBorder = gmx::VectorOperations::ReadAndTransformGeometry(
		oOptionParser.GetOptionValue("-v"), oMercTMS.GetTilingSRSRef());

	
	GeoRasterBuffer* poGeoBuffer = GeoRasterBuffer::InitFromNDVITilesList(listTileContainers, 
		poVecBorder,-1,nZoom);
	
	if (!poGeoBuffer)
	{
		cout << "ERROR: InitFromNDVITilesList\n";
		return 5;
	}

	//debug
	//poGeoBuffer->SaveGeoRefFile("F:\\mpotanin\\byukreevka\\123.tif");
	//end-debug

	int* panIntervals = 0;
	ClassifiedRasterBuffer* poClassifiedBuffer =
		GMXString::MakeLower(oOptionParser.GetOptionValue("-m")) == "intervals" ?
		poGeoBuffer->ClassifyEqualIntervals(nClasses, panIntervals, dblSTDCoeff) :
		poGeoBuffer->ClassifyEqualArea(nClasses, panIntervals);
	

	if (poClassifiedBuffer == 0)
	{
		cout << "ERROR: classify error" << endl;
		return 5;
	}
		
	poClassifiedBuffer->ReplaceByInterpolatedValues(poVecBorder, 1, 1);
	
	for (int i = 0; i < nMajorityFilterNumRun; i++)
		poClassifiedBuffer->ApplyMajorityFilter(nWinSize);
		
	//debug
	

	//poClassifiedBuffer->SaveGeoRefFile("F:\\mpotanin\\byukreevka\\field6_IVI_5class_interpolate_2.tif");
	
	//debug
	///*
	if (oOptionParser.GetOptionValue("-o_geojson") != "")
	{

		map<int, OGRMultiPolygon*> mapZones = poClassifiedBuffer->Polygonize();
		ZoningMap oZM;
		oZM.InitDirectly(mapZones);
		oZM.Clip(poVecBorder);
		
		OGRSpatialReference oWGS84SRS;
		oWGS84SRS.SetWellKnownGeogCS("WGS84");

		oZM.TransformToSRS(oMercTMS.GetTilingSRSRef(), &oWGS84SRS);

		oZM.SaveToVectorFile(oOptionParser.GetOptionValue("-o_geojson"),&oWGS84SRS);

		//oZM.SaveToFile(oOptionParser.GetOptionValue("-o_geojson"));
		//oZM.FilterByArea(5000);
		//oZM.SaveToFile("F:\\mpotanin\\byukreevka\\zones3_clipped_nofilt_field6.shp", "");
	}
	//*/
	//end-debug	


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


	//debug
	poClassifiedBuffer->AdjustExtentToClippedArea();
	//end-debug

	poClassifiedBuffer->set_color_table(poColTab);
	poClassifiedBuffer->SaveGeoRefFile(oOptionParser.GetOptionValue("-o_png") != "" ?
		oOptionParser.GetOptionValue("-o_png") : oOptionParser.GetOptionValue("-o_tif"));





	delete[]panIntervals;
	delete(poClassifiedBuffer);
	delete(poColTab);
	delete(poGeoBuffer);


	//double dblPixelOfset = 2;
	//double dblSTDCoeff = 1.6;
	//int nWinSize = 5;
	//int nClasses

	//
	//Classifier::Init by geometry object
	//static Classifier::ParseVectorFileForBatchProcessing
	//static Classifier::PixelsToPolygons
	//



	/*
	Classifier oClassifier;
	if (!oClassifier.Init(listTileContainers, 
						oOptionParser.GetOptionValue("-v"),
						oOptionParser.GetOptionValue("-z") == "" ? 0 : atoi(oOptionParser.GetOptionValue("-z").c_str())))
	{
		cout << "ERROR: reading input tile containers" << endl;
		return 5;
	}
	

	
	
	GeoRasterBuffer* poClassifiedRaster = 0;

	if (oOptionParser.GetOptionValue("-ranges") == "")
	{
		
	}
	else
	{
		list<string> listRanges = GMXString::SplitCommaSeparatedText(oOptionParser.GetOptionValue("-ranges"));
		int* panRanges = new int[listRanges.size()];
		int i = 0;
		for (string strVal : listRanges)
		{
			panRanges[i] = atoi(strVal.c_str());
			i++;
		}
		poClassifiedRaster = oClassifier.ClassifyByPredefinedIntervals(listRanges.size() - 1, panRanges);
		nClasses = listRanges.size() - 1;
	}
	

	
	if (nMajorityFilterNumRun > 0)
	{
		for (int i = 0; i < nMajorityFilterNumRun; i++)
		{
			Classifier::ApplyMajorityFilter((unsigned char*)poClassifiedRaster->get_pixel_data_ref(),
				poClassifiedRaster->get_x_size(),
				poClassifiedRaster->get_y_size(),
				nClasses,
				nWinSize);
		}
		//poClassifiedRaster->Burn(oOptionParser.GetOptionValue("-v"), 2);
		poClassifiedRaster->Clip(oOptionParser.GetOptionValue("-v"));
	}
	
	


	
	//debug
	poClassifiedRaster->PolygonizeRange(2, 2);

	//end-debug

	if (oOptionParser.GetOptionValue("-o_csv") != "")
	{
		string strOutputText;
		for (string strFile : listTileContainers)
			strOutputText += strFile + "\n";
		GMXFileSys::WriteToTextFile(oOptionParser.GetOptionValue("-o_csv"), strOutputText);
	}
	


	*/

	return 0;
}

