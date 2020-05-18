#pragma once
#include "stdafx.h"
//#include "isoclus_v2/Image.h"

using namespace gmx;

namespace agrigate
{
	class ClassifiedRasterBuffer;


	//ToDo
	//Clip: add buffer offset
	//Polygonize: intersect with clipline
	//
	//
	//

	class GeoRasterBuffer : public gmx::RasterBuffer
	{
	public:
		//PolygonizeRanges
		//VectorOperations::SimplifyPolygons
		//VectorOperations::CensorSmallPolygons

		static GeoRasterBuffer* InitFromNDVIFiles(list<string> listNDVIFiles,
													string strVectorFile = "",
													bool bSaveIVI = true,
													double dblPixelBuffer = 0);
		
		static GeoRasterBuffer* InitFromNDVITiles(list<string> listNDVITiles, 
													OGRGeometry* poVectorMask, 
													bool bSaveIVI=true, 
													bool bMosaicMode = false, 
													double dblPixelBuffer = 0, 
													int nZoom = 0);
		
		//assumed: data_type_==GDT_Int16, single band
		ClassifiedRasterBuffer* ClassifyByPredefinedIntervals(int nNumClass, int* panIntervals);
		
		//assumed: data_type_==GDT_Int16
		ClassifiedRasterBuffer* ClassifyByPredefinedQuantiles(int nNumClass, double* padblQuantiles);

		//assumed: data_type_==GDT_Int16, single band
		ClassifiedRasterBuffer* ClassifyEqualArea(int nNumClass, int* &panIntervals);

		//assumed: data_type_==GDT_Int16, single band
		ClassifiedRasterBuffer* ClassifyEqualIntervals(int nNumClass, int* &panIntervals, double dblSTDCoeff = 0);
	
		
		GeoRasterBuffer(){ m_poSRS = 0; m_nNODV = 0; }
		~GeoRasterBuffer()
		{ 
			if (m_poSRS != 0) OSRRelease(m_poSRS);
		}
		
		bool Clip(OGRGeometry* poClipGeom, double dblPixelBuffer = 0);
		bool Clip(string strVectorFile, double dblPixelOffset = 0);
	
		
		GeoRasterBuffer* BurnVectorMask(OGRGeometry* poClipGeom, double dblPixelOffset = 0);
		
		bool SaveGeoRefFile(string strRasterFile);
		bool CloneGeoRef(GeoRasterBuffer* poSrcBuffer);
		OGRSpatialReference* GetSRSRef() { return m_poSRS; };
		bool SetTMSGeoRef(ITileMatrixSet*  poSRS, int z, int minx, int miny, int maxx, int maxy);
		bool SetGeoRef(OGRSpatialReference* poSRS, OGREnvelope oEnvp, double dblRes);



		GDALDataset* CreateInMemGDALDataset(bool bCopyData = true);
		
		//OGRMultiPolygon* PolygonizeRange(int nMinValue, int nMaxValue);
		
		//assumed: data_type_==GDT_Int16
		map<int, OGRMultiPolygon*> Polygonize(int nOffset, int nStep, int nUpperBound);

		map<int, OGRMultiPolygon*> Polygonize();

		OGREnvelope GetEnvelope(){ return m_oEnvp; };
		double GetPixelSize(){ return m_dblRes; };


		
		//bool CalculateContour(int nMinValue, int nMaxValue, OGRGeometry* &poSRSGeometry);
		//bool CalculateContour(int nMinValue, int nMaxValue, OGRGeometry* &poPixelSpaceGeometry);
		//bool CalculateContours(...);
	protected:


		double* CalcPixelRanks(int nBand = -1);

		
		//template <typename T> GeoRasterBuffer* CreateMaskBufferByRanges(
		//	T type,
		//	int nOffset,
		//	int nStep,
		//	int nUpperBound);

		//template <typename T> GeoRasterBuffer* CreateMaskBufferByRange(T type, int nMinValue, int nMaxValue);
		//OGRGeometry* Clip(OGRGeometry* poClipGeom);
		map<int, OGRMultiPolygon*> Polygonize(GeoRasterBuffer* poGeoBuffer);
		bool CalcGeoTransform(double *padblGeoTransform);
		/*		
			bool GeoRasterBuffer::TraceEdge(int nMinValue,
			int nMaxValue,
			int nX0,
			int nY0,
			int nX1,
			int nY1,
			list<pair<int, int>> &listEdgeLine);
			*/

	protected:
		OGRSpatialReference* m_poSRS;
		double m_dblRes;
		OGREnvelope m_oEnvp;
		int m_nNODV;
	};

	class Classifier
	{
	public:
		//georeference
		Classifier() { m_poGeoBuffer = 0; };
		~Classifier(){ delete(m_poGeoBuffer); };
		
		//bool ClassifyISOCLUSMethod(string strParams,string strBasePath,string strISOFile);
		//input: multiband raster, numclasses //output - classified array
		//init parameters and call method ISOCLUS::Image::RunISOCLUS()

				
		//GeoRasterBuffer* GetGeoRasterBufferRef() { return m_poGeoBuffer; };
		//static bool ApplyMajorityFilter(unsigned char* pabClasses, int w, int h, int nClasses, int nWinSize);
		static bool ConvertPixelsToPolygons(int nZ, string strNDVI, string strOutputVector);
	
	protected:
		GeoRasterBuffer* m_poGeoBuffer;
	};

	
	//assumed: data_type_==GDT_Byte, num_bands_=1
	class ClassifiedRasterBuffer : public GeoRasterBuffer
	{
	public:
		//ToDo
    ClassifiedRasterBuffer* Clone();
		ClassifiedRasterBuffer(int nNumClasses)
		{ 
			m_nClasses = nNumClasses; 
			GeoRasterBuffer();
		};
		
		bool ApplyMajorityFilter(int nWinSize, 
                              bool bOnlyNoDataPixels = false, 
                              unsigned char* pabInterpolationMask = 0);
		int GetNumClasses() { return m_nClasses; };
		bool ReplaceByInterpolatedValues(OGRGeometry* poVector, double dblPixelInward,
			double dblPixelOutward);
		bool AdjustExtentToClippedArea();
		
    //bool PolygonizePixels(string strOutputVectorFile, 
    //                      bool bSaveTo4326 = false);
    bool PolygonizePixels(string strOutputVectorFile, 
                          OGRGeometry* poClipMask, 
                          bool bSaveTo4326 = false);

	protected:
		int m_nClasses;
	};

	class GeoRasterBufferCollection
	{
		OGREnvelope CalculateBundleEnvelope();
		~GeoRasterBufferCollection();
		
		static GeoRasterBufferCollection* InitBundleFromFileInput(string strNDVIFilesTable,
			string strVectorFile, string strBasePath, int nZoom = 13 );
		bool Init(map<string,list<string>> mapNDVITiles, 
			map<string,OGRGeometry*> mapClipGeoms, int nZoom = 13);

		//ToDo
		//Classify
		//PolygonizeToFile

		bool SaveToFile(string strFile, GDALColorTable *poColTab = 0);
	

		map<string, GeoRasterBuffer*> m_mapBuffers;
	};

	struct ImageMeta
	{
		string strFileName;
		int nYear;
		int nDOY;
		double dblMean;
	};


	class NDVIProfile
	{
	public:
		bool ParseInputData(string strHRMean, string strMODISMax, 
							string strAliasPath = "", string strRealPath = "");
		list<string> SelectInputForClassification();
		static bool ParseDateString(string strDate, int* pnYear, int* pnMonth, int* pnDay, int* pnDOY = 0);
		static string GetFullPathBySceneid(string strSceneid, string strBasePath);
	protected:
		map<int, pair<int, double>> m_mapMODISMax;

		double m_padblYearMax[10];
		list<ImageMeta> m_listMetadata;
	};

	class ZoningMap
	{
	public:
		ZoningMap(){};
		ZoningMap(map<int, OGRMultiPolygon*> mapZones)
		{
			InitMakingCopy(mapZones);
		};

		bool TransformToSRS(OGRSpatialReference* poFromSRS, OGRSpatialReference* poToSRS);
		bool InitDirectly(map<int, OGRMultiPolygon*> mapZones);
		bool InitMakingCopy(map<int, OGRMultiPolygon*> mapZones);
		bool Clear();
		~ZoningMap(){ Clear();};

		bool SaveToLayer(OGRLayer* poLayer, string strDNField, 
						string strZoningMapNameField, string strZoningMapName);
		
		bool SaveToVectorFile(string strFileName, OGRSpatialReference* poSRS=0);
		
		bool FilterByArea(double dblThreshold);

		bool Clip(OGRGeometry* poClipGeometry);
	 
	
	protected:
		map<int, OGRPolygon*> FindAllBorderingPolygons(OGRPolygon* poPolygon, double dblThreshold=0);
		bool RemoveEmptyZones();

	protected:
		map<int, OGRMultiPolygon*> m_mapZones;
	};

	class ZoningMapCollection  //ZoningMapCollection
	{
	public:
		ZoningMapCollection(OGRSpatialReference* poSRS)
		{
			m_poSRS = poSRS->Clone();
		};

		~ZoningMapCollection()
		{
			for (auto iter : m_mapZoningCollection)
				delete(iter.second);
			delete(m_poSRS);
		};
		
		bool SaveToVectorFile(string strFileName);
		
		bool AddZoningMap(string strZoningName, ZoningMap* poZM)
		{
			m_mapZoningCollection[strZoningName] = poZM;
			return true;
		};

	protected:
		map<string,ZoningMap*> m_mapZoningCollection;
		OGRSpatialReference* m_poSRS;
	};

	class ClassifiedRastersCollection
	{
	

	};

	class ZoningMapCollectionProcessor
	{
		static map<string, pair<OGRGeometry*, list<string>>>* 
			ParseConsoleInput(string strTextFile, string strVectorFile);
		//...
		
	protected:

	protected:
		map<string,pair<OGRGeometry*,list<string>>> m_mapCollectionData;
	};

	//create CollectionConsoleInput
	//loop geometries and run processing
	//write output vectors into ZoningMapCollection object
	//write output classified rasters into ClassifiedRasterCollection object
}

