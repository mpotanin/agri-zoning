#pragma once
#include "stdafx.h"
#include "isoclus_v2/Image.h"

using namespace gmx;

namespace agrigate
{
	class GeoRasterBuffer : public gmx::RasterBuffer
	{
	public:
		GeoRasterBuffer(){ m_poSRS = 0; }
		~GeoRasterBuffer(){ if (m_poSRS != 0) OSRRelease(m_poSRS); }
		bool SetGeoRef(ITileMatrixSet*  poSRS, int z, int minx, int miny, int maxx, int maxy);
		bool Clip (OGRGeometry* poClipGeom);
		bool Clip(string strVectorFile, double dblPixelOffset = 0);
		bool Polygonize(string strOutputVectorFile);
		
		GeoRasterBuffer* Burn(string strVectorFile, double dblPixelOffset = 0);
		
		bool SaveGeoRefFile(string strOutput);
		GDALDataset* CreateGDALDataset(bool bCopyData = true);
		bool CloneGeoRef(GeoRasterBuffer* poSrcBuffer);
		OGRSpatialReference* GetSRSRef() { return m_poSRS; };
	
	protected:
		bool CalcGeoTransform(double *padblGeoTransform);

	protected:
		OGRSpatialReference* m_poSRS;
		double m_dblRes;
		OGREnvelope m_oEnvp;		
	};

	class Classifier
	{
	public:
		//georeference
		Classifier() { m_poGeoBuffer = 0; };
		~Classifier(){ delete(m_poGeoBuffer); };
		bool Init (list<string> listContainers, string strVectorBorder, int nZoom = 0);
		//GeoRasterBuffer* Classify (params...);
		bool ClassifyISOCLUSMethod(string strParams,string strBasePath,string strISOFile);
		//input: multiband raster, numclasses //output - classified array
		//init parameters and call method ISOCLUS::Image::RunISOCLUS()


		GeoRasterBuffer* ClassifySumMethod(int nNumClass, double dblSTDCoeff, RasterBuffer* poMaskBuffer = 0);
		GeoRasterBuffer* ClassifySumMethod(int nNumClass, int* panRanges);


		GeoRasterBuffer* GetGeoRasterBufferRef() { return m_poGeoBuffer; };
		static bool ApplyMajorityFilter(unsigned char* pabClasses, int w, int h, int nClasses, int nWinSize);
		//static bool DilationFilter(unsigned char* pabClasses, int w, int h, int nClasses, int nWinSize);
	protected:
		uint16_t* CalcIVIRaster();

	protected:
		GeoRasterBuffer* m_poGeoBuffer;
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

}

