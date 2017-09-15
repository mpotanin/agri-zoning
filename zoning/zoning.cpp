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

int nDescriptors = 5;
const GMXOptionDescriptor asDescriptors[] =
{
	{ "-mean", 0, 0, "input mean csv" },
	{ "-max", 0, 0, "input max csv" },
	{ "-v", 0, 0, "input geojson" },
	{ "-o_png", 0, 0, "output png" },
	{ "-o_csv", 0, 0, "output csv file" }
};

int nExamples = 1;
const string astrUsageExamples[] =
{
	"zoning -mean mean.csv -max max.csv -v border.geojson -o_png zones.png -o_csv sceneid_list.csv",
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

	//экстрополяция до границ ?
	//отбор снимков
	//края в MajorituFilter - заполнение пустот
	//сначала cmd

	//string strMeanCSV = "F:\\mpotanin\\FieldStats\\test\\mean.csv";
	//string strMaxCSV = "F:\\mpotanin\\FieldStats\\test\\max.csv";
	//agrigate::NDVIProfile oNP;
	//bool bResult = oNP.ParseInputData(strMeanCSV, strMaxCSV);
	/////////////////////////////////////////////////////////////////////////////
	/*
	string strVectorBorder = "F:\\mpotanin\\FieldStats\\test\\357457.geojson";
	listContainers.push_back("\\\\tinkerbell-smb\\ifs\\kosmosnimki\\Operative\\alt_proc\\ls8\\ndvi\\2015\\2015-06-10\\LC81720282015161LGN00_ndvi.tiles");
	listContainers.push_back("\\\\tinkerbell-smb\\ifs\\kosmosnimki\\Operative\\alt_proc\\ls8\\ndvi\\2016\\2016-04-09\\LC81720292016100LGN00_ndvi.tiles");
	listContainers.push_back("\\\\tinkerbell-smb\\ifs\\kosmosnimki\\Operative\\alt_proc\\ls8\\ndvi\\2017\\2017-04-12\\LC81720282017102LGN00_ndvi.tiles");
	string strOutput = "F:\\mpotanin\\FieldStats\\test\\357457_.png";
	*/

	agrigate::NDVIProfile oNDVIProfile;
	if (!oNDVIProfile.ParseInputData(oOptionParser.GetOptionValue("-mean"),
		oOptionParser.GetOptionValue("-max")))
	{
		cout << "ERROR: can't parse input mean, max csv files" <<endl;
		return 3;
	}

	list<string> listTileContainers = oNDVIProfile.SelectInputForClassification();
	
	if (listTileContainers.size() == 0)
	{
		cout << "ERROR: no image fits criterias" << endl;
		return 4;
	}

	
	int nWinSize = 5;
	int nClasses = 5;
	int nRepeatFilter = 2;
	GDALColorTable *poColTab = new GDALColorTable(GPI_RGB);
	GDALColorEntry arrColors[6];
	arrColors[0].c1 = 0; arrColors[0].c2 = 0; arrColors[0].c3 = 0;
	arrColors[1].c1 = 255; arrColors[1].c2 = 0; arrColors[1].c3 = 0;
	arrColors[2].c1 = 247; arrColors[2].c2 = 209; arrColors[2].c3 = 59;
	arrColors[3].c1 = 212; arrColors[3].c2 = 255; arrColors[3].c3 = 190;
	arrColors[4].c1 = 76; arrColors[4].c2 = 227; arrColors[4].c3 = 0;
	arrColors[5].c1 = 47; arrColors[5].c2 = 140; arrColors[5].c3 = 30;
	for (int i = 0; i <= nClasses; i++)
		poColTab->SetColorEntry(i, &arrColors[i]);
	/////////////////////////////////////////////////////////////////////////////

	Classifier oClassifier;
	if (!oClassifier.Init(listTileContainers, oOptionParser.GetOptionValue("-v")))
	{
		cout << "ERROR: reading input tile containers" << endl;
		return 5;
	}

	GeoRasterBuffer* poGeoBuffer = oClassifier.ClassifySumMethod();
	for (int i = 0; i < 2; i++)
	{
		Classifier::ApplyMajorityFilter((unsigned char*)poGeoBuffer->get_pixel_data_ref(),
			poGeoBuffer->get_x_size(),
			poGeoBuffer->get_y_size(),
			nClasses,
			nWinSize);
	}
	
	poGeoBuffer->set_color_table(poColTab);
	poGeoBuffer->SaveGeoRefFile(oOptionParser.GetOptionValue("-o_png"));

	string strOutputText;
	for (string strFile : listTileContainers)
		strOutputText += strFile + "\n";
	GMXFileSys::WriteToTextFile(oOptionParser.GetOptionValue("-o_csv"), strOutputText);

	//png 32bit

	delete(poGeoBuffer);
	delete(poColTab);


	return 0;
}

