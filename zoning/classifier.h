#pragma once
#include "stdafx.h"

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
		bool SaveGeoRefFile(string strOutput);
		GDALDataset* CreateGDALDataset(bool bCopyData = true);
		bool CloneGeoRef(GeoRasterBuffer* poBuffer);
	
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
		bool Init (list<string> listContainers, string strVectorBorder);
		//GeoRasterBuffer* Classify (params...);
		GeoRasterBuffer* ClassifySumMethod(int nNumClass = 5, double dblSTDCoeff = 1.6);
		static bool ApplyMajorityFilter(unsigned char* pabClasses, int w, int h, int nClasses, int nWinSize);
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
		bool ParseInputData(string strHRMean, string strMODISMax);
		list<string> SelectInputForClassification();
		static bool ParseDateString(string strDate, int* pnYear, int* pnMonth, int* pnDay, int* pnDOY = 0);
	protected:
		map<int, pair<int, double>> m_mapMODISMax;

		double m_padblYearMax[10];
		list<ImageMeta> m_listMetadata;
	};

}

