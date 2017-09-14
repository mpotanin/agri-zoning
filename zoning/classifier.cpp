#include "stdafx.h"
#include "classifier.h"


namespace agrigate
{
	bool GeoRasterBuffer::CloneGeoRef(GeoRasterBuffer* poBuffer)
	{
		if (!poBuffer) return false;
		if (!poBuffer->m_poSRS) return false;
		m_poSRS = poBuffer->m_poSRS->Clone();
		
		m_oEnvp = poBuffer->m_oEnvp;
		m_dblRes = poBuffer->m_dblRes;

		return true;
	}

	bool GeoRasterBuffer::SetGeoRef(ITileMatrixSet*  poSRS, int z, int minx, int miny, int maxx, int maxy)
	{
		m_poSRS = poSRS->GetTilingSRSRef()->Clone();
		m_dblRes = poSRS->CalcPixelSizeByZoom(z);
		m_oEnvp = poSRS->CalcEnvelopeByTileRange(z, minx, miny, maxx, maxy);
		
		return false;
	}

	GDALDataset* GeoRasterBuffer::CreateGDALDataset(bool bCopyData)
	{
		string strTiffInMem = "/vsimem/tiffinmem_" + GMXString::ConvertIntToString(rand());
		double padblGeoTransform[6];
		if (!CalcGeoTransform(padblGeoTransform)) return 0;
		
		GDALDataset*	poVrtDS = (GDALDataset*)GDALCreate(
			GDALGetDriverByName("GTiff"),
			strTiffInMem.c_str(),
			x_size_,
			y_size_,
			this->num_bands_,
			this->data_type_,
			0
			);

		if (bCopyData) 
			poVrtDS->RasterIO(GF_Write, 0, 0, x_size_, y_size_, get_pixel_data_ref(),
							x_size_, y_size_, data_type_,num_bands_, 0, 0, 0, 0);

		poVrtDS->SetGeoTransform(padblGeoTransform);
		char *pachWKT;
		m_poSRS->exportToWkt(&pachWKT);
		poVrtDS->SetProjection(pachWKT);
		OGRFree(pachWKT);
		
		return poVrtDS;
	}

	bool GeoRasterBuffer::CalcGeoTransform(double *padblGeoTransform)
	{
		if (m_poSRS == 0 || padblGeoTransform == 0) return false;

		padblGeoTransform[5] = -(padblGeoTransform[1] = m_dblRes);
		padblGeoTransform[0] = m_oEnvp.MinX;
		padblGeoTransform[3] = m_oEnvp.MaxY;
		padblGeoTransform[2] = (padblGeoTransform[4] = 0);

		return true;
	}


	bool GeoRasterBuffer::Clip(OGRGeometry* poClipGeom)
	{
		double padblGeoTransform[6];
		if (!CalcGeoTransform(padblGeoTransform)) return false;
				
		GDALDataset* poBufferDS = CreateGDALDataset();
		GDALDataset* poClippedDS = CreateGDALDataset(false);

		GDALWarpOptions *poWarpOptions = GDALCreateWarpOptions();
		poWarpOptions->papszWarpOptions = 0;
		poWarpOptions->hSrcDS = poBufferDS;
		poWarpOptions->hDstDS = poClippedDS;
		poWarpOptions->dfWarpMemoryLimit = 250000000;
		double			dblError = 0.125;
		
		char* pachSRSWKT;
		m_poSRS->exportToWkt(&pachSRSWKT);
		poWarpOptions->hCutline = VectorOperations::ConvertFromSRSToPixelLine(poClipGeom, padblGeoTransform);
		poWarpOptions->pTransformerArg = GDALCreateApproxTransformer(GDALGenImgProjTransform,
			GDALCreateGenImgProjTransformer(poBufferDS, pachSRSWKT, poClippedDS, pachSRSWKT, false, 0.0, 1),
			dblError);
		poWarpOptions->pfnTransformer = GDALApproxTransform;
		//poWarpOptions->eResampleAlg = GRA_NearestNeighbour;
		poWarpOptions->eResampleAlg = GRA_Cubic;

		// Initialize and execute the warp operation. 
		GDALWarpOperation gdal_warp_operation;
		gdal_warp_operation.Initialize(poWarpOptions);

		bool  warp_error = (CE_None != gdal_warp_operation.ChunkAndWarpImage(0, 0, x_size_, y_size_));

		GDALDestroyApproxTransformer(poWarpOptions->pTransformerArg);
		delete((OGRGeometry*)poWarpOptions->hCutline);
		poWarpOptions->hCutline = 0;
		GDALDestroyWarpOptions(poWarpOptions);
		OGRFree(pachSRSWKT);
		poClippedDS->RasterIO(GF_Read, 0, 0, x_size_, y_size_, get_pixel_data_ref(),
			x_size_, y_size_, data_type_, num_bands_, 0, 0, 0, 0);
		GDALClose(poBufferDS);
		GDALClose(poClippedDS);
				

		return true;
	}

	bool GeoRasterBuffer::SaveGeoRefFile(string strOutput)
	{
		SaveBufferToFile(strOutput);
		GMXFileSys::WriteWLDFile(strOutput, m_oEnvp.MinX, m_oEnvp.MaxY, m_dblRes);
		
		return true;
	}

	bool Classifier::ApplyMajorityFilter(unsigned char* pabClasses, int w, int h, int nClasses, int nWinSize)
	{
		int n, k;
		int *panFreq = new int[nClasses + 1];
		int nWinSize_2 = nWinSize / 2;
		int nMaxFreqClass;

		unsigned char* pabClassesFiltered = new unsigned char[w*h];
		for (int y = 0; y<h; y++)
		for (int x = 0; x<w; x++)
			pabClassesFiltered[w*y + x] = pabClasses[y*w + x];


		for (int y = nWinSize_2; y < (h - nWinSize_2); y++)
		{
			for (int x = nWinSize_2; x < (w - nWinSize_2); x++)
			{
				for (int q = 0; q <= nClasses; q++)
					panFreq[q] = 0;

				for (int i = -nWinSize_2; i <= nWinSize_2; i++)
				{
					for (int j = -nWinSize_2; j <= nWinSize_2; j++)
					{
						n = w*(y + i) + x + j;
						panFreq[pabClasses[n]]++;
					}
				}

				nMaxFreqClass = pabClasses[w*y + x];
				for (int q = 0; q <= nClasses; q++)
					nMaxFreqClass = panFreq[q] > panFreq[nMaxFreqClass] ? q : nMaxFreqClass;
				pabClassesFiltered[w*y + x] = nMaxFreqClass;
			}
		}

		for (int y = 0; y<h; y++)
		for (int x = 0; x<w; x++)
			pabClasses[y*w + x] = pabClassesFiltered[w*y + x];

		delete[]panFreq;
		delete[]pabClassesFiltered;

		return true;
	}


	GeoRasterBuffer* Classifier::ClassifySumMethod(int nNumClass, double dblSTDCoeff)
	{
		int nWidth = m_poGeoBuffer->get_x_size();
		int nHeight = m_poGeoBuffer->get_y_size();
		uint16_t* panIVI = new uint16_t[nWidth * nHeight];
		
		int nBands = m_poGeoBuffer->get_num_bands();
		unsigned char* panNDVI = (unsigned char*)m_poGeoBuffer->get_pixel_data_ref();
		int t;
		int n = nWidth*nHeight;
		for (int i = 0; i < nHeight; i++)
		{
			for (int j = 0; j < nWidth; j++)
			{
				t = i*nWidth + j;
				panIVI[t] = 0;
				for (int b = 0; b < nBands; b++)
					panIVI[t] += panNDVI[b*n + t];
			}
		}

		
		double dblMean = 0;
		double dblSTD = 0;
		int nNumPixels = 0;

		for (int i = 0; i < nHeight; i++)
		{
			for (int j = 0; j < nWidth; j++)
			{
				t = i*nWidth + j;
				if (panIVI[t] != 0)
				{
					nNumPixels++;
					dblMean += panIVI[t];
					dblSTD += panIVI[t] * panIVI[t];
				}
			}
		}

		dblMean /= nNumPixels;
		dblSTD = sqrt((dblSTD / nNumPixels) - dblMean*dblMean);
		double dblMin = dblMean - dblSTD * dblSTDCoeff;
		double dblMax = dblMean + dblSTD * dblSTDCoeff;
		double dblDif = (dblMax - dblMin) / nNumClass;

		unsigned char* panClasses = new unsigned char[n];
		for (int i = 0; i < nHeight; i++)
		{
			for (int j = 0; j < nWidth; j++)
			{
				t = i*nWidth + j;
				panClasses[t] = (panIVI[t] == 0) ? 0 :
					panIVI[t] > dblMax ? nNumClass :
					panIVI[t] < dblMin ? 1 :
					(int)((panIVI[t] - dblMin - 1e-5) / dblDif) + 1;
			}	
		}
		GeoRasterBuffer *poOutputBuffer = new GeoRasterBuffer();
		poOutputBuffer->CreateBuffer(1, nWidth, nHeight, panClasses);
		
		poOutputBuffer->CloneGeoRef(m_poGeoBuffer);
		

		delete[]panClasses;
		delete[]panIVI;

		return poOutputBuffer;
	}


	bool Classifier::Init(list<string> listContainers, string strVectorBorder)
	{
		if (listContainers.size() == 0) return false;

		ITileContainer* poContainer =  TileContainerFactory::OpenForReading(*listContainers.begin());
		if (poContainer == 0) return 0;
		MercatorTileMatrixSet oTMS(poContainer->GetProjType());
		int z = poContainer->GetMaxZoom();
		poContainer->Close();
		delete(poContainer);
		
		OGRGeometry* poFieldGeometry = VectorOperations::ReadAndTransformGeometry(strVectorBorder,
																				oTMS.GetTilingSRSRef());
		OGREnvelope oFieldEnvp;
		poFieldGeometry->getEnvelope(&oFieldEnvp);
		
		int minx, maxx, miny, maxy;
		oTMS.CalcTileRange(oFieldEnvp, z, minx, miny, maxx, maxy);

		m_poGeoBuffer = new GeoRasterBuffer();
		m_poGeoBuffer->CreateBuffer(listContainers.size(), 256 * (maxx - minx + 1), 256 * (maxy - miny + 1));
		
		int f = -1;
		for (string strContainer : listContainers)
		{
			f++;
			poContainer = TileContainerFactory::OpenForReading(strContainer);
			if (poContainer == 0)
			{
				delete (poFieldGeometry);
				return false;
			}

			for (int x = minx; x <= maxx; x++)
			{
				for (int y = miny; y <= maxy; y++)
				{
					unsigned char* pabyTileData = 0;
					unsigned int nSize = 0;

					if (!poContainer->GetTile(z, x, y, pabyTileData, nSize))
					{
						delete (poFieldGeometry);
						return false;
					}

					RasterBuffer oTileBuffer;
					oTileBuffer.CreateBufferFromPngData(pabyTileData, nSize);
					m_poGeoBuffer->SetPixelDataBlock(256*(x-minx), 256*(y-miny), 256, 256, oTileBuffer.get_pixel_data_ref(), f, f);
					delete[]pabyTileData;
				}
			}
			poContainer->Close();
			delete(poContainer);
		}
		
		m_poGeoBuffer->SetGeoRef(&oTMS, z, minx, miny, maxx, maxy);
		m_poGeoBuffer->Clip(poFieldGeometry);
		delete(poFieldGeometry);

		return true;
	}

	list<string>  NDVIProfile::SelectInputForClassification()
	{
		list<string> listBestImages;

		//loop ordered list
		//algorithm1
		//for each year max_value >= 0.5
		//select first image
		//if date earlier than max than select
		int nLastYear = 0;
		for (ImageMeta sMeta : m_listMetadata)
		{
			if (sMeta.nYear == nLastYear) continue;
			else if (m_mapMODISMax[sMeta.nYear].second < 0.5) continue;
			else if (sMeta.nDOY>m_mapMODISMax[sMeta.nYear].first) continue;
			else
			{
				nLastYear = sMeta.nYear;
				listBestImages.push_back(sMeta.strFileName);
			}
		}

		return listBestImages;
	}

	bool ParseDateString(string strDate, int* pnYear, int* pnMonth, int* pnDay, int* pnDOY)
	{
		regex regYYYYMMDD("\\d{7}");
		regex regYYYY_MM_DD("\\d{4}\\.\\d{2}\\.\\d{2}");
		regex regDD_MM_YYYY("\\d{2}\\.\\d{2}\\.\\d{4}");
		
		int panLeapYear[12] = { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
		int panSimpleYear[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };

		string strYYYY;
		string strMM;
		string strDD;
		if (regex_match(strDate, regYYYYMMDD))
		{
			strYYYY = strDate.substr(0, 4);
			strMM = strDate.substr(4, 2);
			strDD = strDate.substr(6, 2);
		}
		else if (regex_match(strDate, regYYYY_MM_DD))
		{
			strYYYY = strDate.substr(0, 4);
			strMM = strDate.substr(5, 2);
			strDD = strDate.substr(8, 2);
		}
		else if (regex_match(strDate, regDD_MM_YYYY))
		{
			strYYYY = strDate.substr(6, 4);
			strMM = strDate.substr(3, 2);
			strDD = strDate.substr(0, 2);
		}
		else return false;

		(*pnYear) = atoi(strYYYY.c_str());
		(*pnMonth) = atoi(strMM.c_str());
		(*pnDay) = atoi(strDD.c_str());
		(*pnDOY) = (*pnYear) % 4 == 0 ? panLeapYear[(*pnMonth) - 1] + (*pnDay) :
										panSimpleYear[(*pnMonth) - 1] + (*pnDay);
		return true;
	}


	bool NDVIProfile::ParseInputData(string strHRMeanFile, string strMODISMaxFile)
	{
		list<string> listLines;
		if (!GMXFileSys::ReadTextFile(strHRMeanFile, listLines)) return false;
		listLines.pop_front();
		int nDay, nMonth, nYear, nDOY;
		for (string strLine : listLines)
		{
			if (!NDVIProfile::ParseDateString(strLine.substr(0,10),&nYear,&nMonth,&nDay,&nDOY))
			{
				m_listMetadata.clear();
				return false;
			}
			ImageMeta oMeta;
			oMeta.nYear = nYear;
			oMeta.nDOY = nDOY;
			
			string strFilePath = strLine.substr(19, (strLine.find(';', 19) - 19));
			oMeta.strFileName = strFilePath.find("\\\\tinkerbell-smb\\ifs\\") == string::npos ?
				("\\\\tinkerbell-smb\\ifs\\" + strFilePath) : strFilePath;
			oMeta.dblMean = atof(strLine.substr(strLine.find(';', 19)+1).replace(1,1,".").c_str());
			m_listMetadata.push_back(oMeta);
		}

			
		listLines.clear();
		
		if (!GMXFileSys::ReadTextFile(strMODISMaxFile, listLines))
		{
			m_listMetadata.clear();
			return false;
		}					
		
		listLines.pop_front();
		for (string strLine : listLines)
		{
			if (strLine.size() < 10) break;
			else strLine += ";";

			int nYear = atoi(strLine.substr(0,strLine.find(";")).c_str());
			strLine = strLine.substr(strLine.find(";") + 1);
			double dblMax = atof(strLine.substr(0, strLine.find(";")).c_str());
			strLine = strLine.substr(strLine.find(";") + 1);
			if (!NDVIProfile::ParseDateString(strLine.substr(0, strLine.find(";")), &nYear, &nMonth, &nDay, &nDOY))
			{
				m_listMetadata.clear();
				m_mapMODISMax.clear();
			}
			m_mapMODISMax[nYear] = pair<int,double>(nDOY,dblMax);
		}

		return true;
	}
}



