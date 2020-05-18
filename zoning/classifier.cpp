#include "stdafx.h"
#include "classifier.h"


namespace agrigate
{
	GeoRasterBuffer* GeoRasterBuffer::InitFromNDVIFiles(list<string> listNDVIFiles,
														string strVectorFile,
														bool bSaveIVI,
														double dblPixelBuffer)
	{
		if (listNDVIFiles.size() == 0) return 0;

		string strFirstFile = (*listNDVIFiles.begin());
		gmx::RasterFile oRF;
		if (!oRF.Init(strFirstFile))
		{
			cout << "ERROR: reading file: " << strFirstFile << endl;
			return 0;
		}

		OGRSpatialReference oBufferSRS;
		if (!oRF.GetSRS(oBufferSRS))
		{
			cout << "ERROR: reading spatial reference from file: " << strFirstFile << endl;
			return 0;

		}

		double padblGeoTransform[6];
		oRF.GetGeoTransform(padblGeoTransform);
		
		double dblPixelRes = padblGeoTransform[1];
		
		int nWidth, nHeight;
		oRF.GetPixelSize(nWidth, nHeight);

		OGREnvelope oEnvp;
		oEnvp.MinX = padblGeoTransform[0];
		oEnvp.MaxY = padblGeoTransform[3];
		oEnvp.MaxX = oEnvp.MinX + dblPixelRes * nWidth;
		oEnvp.MinY = oEnvp.MaxY - dblPixelRes * nHeight;

		int nOffsetX = 0;
		int nOffsetY = 0;

		OGRGeometry* poMultiPoly = 0;

		if (strVectorFile != "")
		{
			poMultiPoly = gmx::VectorOperations::ReadIntoSingleMultiPolygon(strVectorFile, &oBufferSRS);
			if (!poMultiPoly)
			{
				cout << "ERROR: reading geometry from file: " << strVectorFile << endl;
				return 0;
			}

			poMultiPoly->getEnvelope(&oEnvp);//TODO: pixelbuffer
			nWidth = ceil(((oEnvp.MaxX - oEnvp.MinX) / dblPixelRes) + 2);
			nHeight = ceil(((oEnvp.MaxY - oEnvp.MinY) / dblPixelRes) + 2);


			oEnvp.MinX = max(padblGeoTransform[0],
				padblGeoTransform[0] + dblPixelRes * (floor((oEnvp.MinX - padblGeoTransform[0]) / dblPixelRes) - 1));
			oEnvp.MaxY = min(padblGeoTransform[3],
				padblGeoTransform[3] - dblPixelRes * (floor((padblGeoTransform[3] - oEnvp.MaxY) / dblPixelRes) - 1));

			oEnvp.MaxX = oEnvp.MinX + nWidth * dblPixelRes;
			oEnvp.MinY = oEnvp.MaxY - nHeight * dblPixelRes;


		}
		
		GeoRasterBuffer *poOutputBuffer = new GeoRasterBuffer();
		bSaveIVI = false; //TODO: if bSaveIVI == true

		poOutputBuffer->CreateBuffer(listNDVIFiles.size(), nWidth, nHeight, 0, GDT_Int16);
		poOutputBuffer->SetGeoRef(&oBufferSRS, oEnvp, dblPixelRes);
		GDALDataset* poVrtDS = poOutputBuffer->CreateInMemGDALDataset(false);
		
		int nCount = 1;

		for (auto strFile : listNDVIFiles)
		{
			//warp from strFile into listNDVIFiles
			GDALDataset* poInputFileDS = (GDALDataset*)GDALOpen(strFile.c_str(), GA_ReadOnly);
			GDALWarpOptions *poWarpOptions = GDALCreateWarpOptions();
			poWarpOptions->papszWarpOptions = 0;

			poWarpOptions->hSrcDS = poInputFileDS;
			poWarpOptions->hDstDS = poVrtDS;

			poWarpOptions->dfWarpMemoryLimit = 250000000;

			double dblErr = 0.125;

			poWarpOptions->nBandCount = 1;
			poWarpOptions->panSrcBands = new int[1];
			poWarpOptions->panDstBands = new int[1];
			poWarpOptions->panSrcBands[0] = 1;
			poWarpOptions->panDstBands[0] = nCount;

			poWarpOptions->padfSrcNoDataReal = new double[1];
			poWarpOptions->padfSrcNoDataImag = new double[1];
			poWarpOptions->padfSrcNoDataReal[0] = poOutputBuffer->m_nNODV;
			poWarpOptions->padfSrcNoDataImag[0] = 0;
			
			poWarpOptions->pfnProgress = GMXPrintProgressStub;
		
			char* pachBufferWKT;
			oBufferSRS.exportToWkt(&pachBufferWKT);

			//char* pachFileWKT;
			//poFileDS->GetProjectionRef()

			poWarpOptions->pTransformerArg = GDALCreateApproxTransformer(GDALGenImgProjTransform,
										GDALCreateGenImgProjTransformer(poInputFileDS, 0, poVrtDS, pachBufferWKT, false, 0.0, 1),
										dblErr);
			poWarpOptions->pfnTransformer = GDALApproxTransform;
			poWarpOptions->eResampleAlg = GRA_Cubic;


			/*
			poWarpOptions->padfSrcNoDataReal = new double[num_bands_];
			poWarpOptions->padfSrcNoDataImag = new double[num_bands_];
			for (int i = 0; i < num_bands_; i++)
			{
				poWarpOptions->padfSrcNoDataReal[i] = m_nNODV;
				poWarpOptions->padfSrcNoDataImag[i] = 0;
			}
			*/

			// Initialize and execute the warp operation. 
			GDALWarpOperation gdal_warp_operation;
			gdal_warp_operation.Initialize(poWarpOptions);

			bool  warp_error = (CE_None != gdal_warp_operation.ChunkAndWarpImage(0, 0, nWidth, nHeight));

			delete[]poWarpOptions->panSrcBands;
			delete[]poWarpOptions->panDstBands;
			poWarpOptions->panSrcBands = 0;
			poWarpOptions->panDstBands = 0;
			delete[]poWarpOptions->padfSrcNoDataReal;
			delete[]poWarpOptions->padfSrcNoDataImag;
			poWarpOptions->padfSrcNoDataReal = 0;
			poWarpOptions->padfSrcNoDataImag = 0;

			GDALDestroyApproxTransformer(poWarpOptions->pTransformerArg);
			GDALDestroyWarpOptions(poWarpOptions);
			OGRFree(pachBufferWKT);
			GDALClose(poInputFileDS);
			
			nCount++;
		}


		poVrtDS->RasterIO(GF_Read, 0, 0, nWidth, nHeight, poOutputBuffer->get_pixel_data_ref(),
			nWidth, nHeight, poOutputBuffer->get_data_type(), poOutputBuffer->get_num_bands(), 0, 0, 0, 0);
		
		if (poMultiPoly)
		{
			poOutputBuffer->Clip(poMultiPoly, dblPixelBuffer);
			delete(poMultiPoly);
		}

		return poOutputBuffer;
		
	}

	bool GeoRasterBuffer::SetGeoRef(OGRSpatialReference* poSRS, OGREnvelope oEnvp, double dblRes)
	{
		if (!poSRS) return false;
		
		m_poSRS = poSRS->Clone();
		m_dblRes = dblRes;
		m_oEnvp = oEnvp;

		return true;
	}

	GeoRasterBuffer*  GeoRasterBuffer::InitFromNDVITiles(
										list<string> listNDVITiles, 
										OGRGeometry* poVectorMask,
										bool bSaveIVI,
										bool bMosaicMode,
										double dblPixelBuffer,
										int nZoom )
	{		
		if (listNDVITiles.size() == 0) return false;
		if (!poVectorMask) return false;

		int nNODV = -32767;
		int* panSign = new int[listNDVITiles.size()];
		int ind = -1;
		for (list<string>::iterator iter = listNDVITiles.begin(); iter != listNDVITiles.end(); iter++)
		{
			ind++;
			if ((*iter)[0] == '-')
			{
				//nNODV = -32767;
				(*iter) = (*iter).substr(1);
				panSign[ind] = -1;
			}
			else panSign[ind] = 1;
		}
		
		GMXTileContainer* poContainer = 
			(GMXTileContainer*)TileContainerFactory::OpenForReading(*listNDVITiles.begin());
		if (poContainer == 0)
		{
			delete[]panSign;
			return 0;
		}
		MercatorTileMatrixSet oWebMercTMS(WEB_MERCATOR);
		
		nZoom = (nZoom == 0) ? poContainer->GetMaxZoom() : nZoom;
		poContainer->Close();
		delete(poContainer);

		OGREnvelope oVectorMaskEnvp;
		poVectorMask->getEnvelope(&oVectorMaskEnvp);
		double dblEnvpOffset = 2 * oWebMercTMS.CalcPixelSizeByZoom(nZoom);
		oVectorMaskEnvp.MinX -= dblEnvpOffset;
		oVectorMaskEnvp.MinY -= dblEnvpOffset;
		oVectorMaskEnvp.MaxX += dblEnvpOffset;
		oVectorMaskEnvp.MaxY += dblEnvpOffset;
		
		int minx, maxx, miny, maxy;
		oWebMercTMS.CalcTileRange(oVectorMaskEnvp, nZoom, minx, miny, maxx, maxy);

		int nWidth = 256 * (maxx - minx + 1);
		int nHeight = 256 * (maxy - miny + 1);
		GeoRasterBuffer *poOutputBuffer =  new GeoRasterBuffer();
		
		bSaveIVI = bMosaicMode ? true : bSaveIVI;

		if (bSaveIVI) poOutputBuffer->CreateBuffer(1, nWidth, nHeight, 0, GDT_Int16);
		else poOutputBuffer->CreateBuffer(listNDVITiles.size(), nWidth, nHeight, 0, GDT_Int16);
		int nBands = poOutputBuffer->get_num_bands();

		poOutputBuffer->InitByValue(poOutputBuffer->m_nNODV = nNODV);
		
		int16_t* panOutputVals = (int16_t*)poOutputBuffer->get_pixel_data_ref();

		ind = -1;
		for (string strIter : listNDVITiles)
		{
			ind++;
			
			if (!(poContainer = (GMXTileContainer*)TileContainerFactory::OpenForReading(strIter)))
			{
				delete(poOutputBuffer);
				cout << "ERROR: reading file " << strIter << endl;
				return false;
			}

			for (int x = minx; x <= maxx; x++)
			{
				for (int y = miny; y <= maxy; y++)
				{
					unsigned char* pabyTileData = 0;
					unsigned int nSize = 0;

					if (!poContainer->GetTile(nZoom, x, y, pabyTileData, nSize))
					{
						if (!bMosaicMode)
						{
						  printf("ERROR: reading tile (%d, %d, %d) from file %s\n", nZoom, x, y, strIter);
						  poContainer->Close();
						  delete(poOutputBuffer);
						  return false;
						}
						else continue;
					}
    
					RasterBuffer oTileBuffer;
					oTileBuffer.CreateBufferFromInMemoryData(pabyTileData, nSize, PNG_TILE);
					delete[]pabyTileData;

					unsigned char* panNDVIVals = (unsigned char*)oTileBuffer.get_pixel_data_ref();
					int nOffsetTop = 256 * (y - miny);
					int nOffsetLeft = 256 * (x - minx);

					int n1, n2;
					for (int i = 0; i < 256; i++)
					{
						for (int j = 0; j < 256; j++)
						{
							n1 = 256 * j + i;
							if (panNDVIVals[n1] != 0)
							{
								n2 = bSaveIVI ? nWidth*(nOffsetTop + j) + nOffsetLeft + i :
												ind*nWidth*nHeight + nWidth * (nOffsetTop + j) + nOffsetLeft + i;
 								if (!bMosaicMode)
								{
									panOutputVals[n2] = panOutputVals[n2] == nNODV ? 0 : panOutputVals[n2];
									panOutputVals[n2] += panSign[ind] > 0 ? (panNDVIVals[n1] - 101) :
															-(panNDVIVals[n1] - 101);
								}
								else
								{
									panOutputVals[n2] = panOutputVals[n2] < (panNDVIVals[n1]-101) ? 
										panNDVIVals[n1]-101 : panOutputVals[n2];
								}
							}
						}
					}
				
				}
			}
			poContainer->Close();
			delete(poContainer);
		}

		poOutputBuffer->SetTMSGeoRef(&oWebMercTMS, nZoom, minx, miny, maxx, maxy);
		poOutputBuffer->Clip(poVectorMask,dblPixelBuffer);
		
		delete[]panSign;
		return poOutputBuffer;
	}

	bool GeoRasterBuffer::CloneGeoRef(GeoRasterBuffer* poBuffer)
	{
		if (!poBuffer) return false;
		if (!poBuffer->m_poSRS) return false;
		m_poSRS = poBuffer->m_poSRS->Clone();
		
		m_oEnvp = poBuffer->m_oEnvp;
		m_dblRes = poBuffer->m_dblRes;

		return true;
	}

	bool GeoRasterBuffer::SetTMSGeoRef(ITileMatrixSet*  poSRS, int z, int minx, int miny, int maxx, int maxy)
	{
		m_poSRS = poSRS->GetTilingSRSRef()->Clone();
		m_dblRes = poSRS->CalcPixelSizeByZoom(z);
		m_oEnvp = poSRS->CalcEnvelopeByTileRange(z, minx, miny, maxx, maxy);
		
		return false;
	}

	/*
	template <typename T> 
	GeoRasterBuffer* GeoRasterBuffer::CreateMaskBufferByRanges(T type, int nOffset, int nStep, int nUpperBound)
	{
		int n = x_size_*y_size_;
		unsigned char* panMaskPixels = new unsigned char[n];
		T* panPixels = (T*)p_pixel_data_;
		for (int i = 0; i < n; i++)
		{
			if ((panPixels[i] < nOffset) || (panPixels[i]>nUpperBound))
				panMaskPixels[i] = m_nNODV;
			else
				panMaskPixels[i] = 1 + ((panPixels[i]-nOffset)/nStep);
		}

		GeoRasterBuffer* poMaskBuffer = new GeoRasterBuffer();
		poMaskBuffer->CreateBuffer(1, x_size_, y_size_, panMaskPixels, GDT_Byte);
		poMaskBuffer->CloneGeoRef(this);
		delete[]panMaskPixels;
		return poMaskBuffer;
	}

	
	template <typename T>
	GeoRasterBuffer* GeoRasterBuffer::CreateMaskBufferByRange(T type, int nMinValue, int nMaxValue)
	{
		int n = x_size_*y_size_;
		int16_t* panMaskPixels = new int16_t[n];
		
				
		T* panPixels = (T*)p_pixel_data_;
		for (int i = 0; i < n; i++)
		{
			panMaskPixels[i] = (panPixels[i] >= nMinValue &&  panPixels[i] <= nMaxValue) ?
				panPixels[i] : 0;
		}
		
		GeoRasterBuffer* poMaskBuffer = new GeoRasterBuffer();
		poMaskBuffer->CreateBuffer(1, x_size_, y_size_, panMaskPixels, GDT_Int16);
		poMaskBuffer->CloneGeoRef(this);
		delete[]panMaskPixels;
		return poMaskBuffer;
	}
	*/

	GDALDataset* GeoRasterBuffer::CreateInMemGDALDataset(bool bCopyData)
	{
		string strTiffInMem = "/vsimem/tiffinmem_" + GMXString::ConvertIntToString(rand());
		double padblGeoTransform[6];
		if (!CalcGeoTransform(padblGeoTransform)) return 0;
		
		GDALDataset*	poVrtDS = (GDALDataset*)GDALCreate(
			GDALGetDriverByName("GTiff"),
			strTiffInMem.c_str(),
			x_size_,
			y_size_,
			num_bands_,
			data_type_,
			0
			);

		for (int b = 1; b <= num_bands_; b++)
			poVrtDS->GetRasterBand(b)->SetNoDataValue(m_nNODV);


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

	bool GeoRasterBuffer::Clip(string strVectorFile, double dblPixelOffset)
	{
		OGRGeometry* poClipGeometry =
      VectorOperations::ReadIntoSingleMultiPolygon(strVectorFile, GetSRSRef());
		if (!poClipGeometry) return false;

		bool bResult = Clip(poClipGeometry,dblPixelOffset);

		delete(poClipGeometry);

		return bResult;
	}

  
	bool  ClassifiedRasterBuffer::AdjustExtentToClippedArea()
	{
		int nTop=0, nLeft=0, nRight=0, nBottom=0;

		//calculate pixel offsets: top, bottom, left, right
		unsigned char* panPixels = (unsigned char*)p_pixel_data_;
		int nArea = x_size_*y_size_;
		int i,j;
		for (i = 0; i<nArea; i++)
			if (panPixels[i] != 0) break;
		nTop = i / x_size_;

		for (i = nArea-1; i>-1; i--)
			if (panPixels[i] != 0) break;
		nBottom = y_size_ -1 - (i / x_size_);

		for (i = 0; i < x_size_; i++)
		{
			for (j = 0; j < y_size_; j++)
			{
				if (panPixels[j*x_size_ + i] != 0) break;
			}
			if (j < y_size_)
			{
				nLeft = i;
				break;
			}
		}

		for (i = x_size_ - 1; i > -1; i--)
		{
			for (j = 0; j < y_size_; j++)
			{
				if (panPixels[j*x_size_ + i] != 0) break;
			}
			if (j < y_size_)
			{
				nRight = x_size_ - 1 - i;
				break;
			}
		}
		


		void* panPixelBlock = 
			GetPixelDataBlock(nLeft, nTop, x_size_ - nLeft - nRight, y_size_ - nTop - nBottom);

		ClearBuffer();

		p_pixel_data_ = panPixelBlock;

		x_size_ -= (nLeft + nRight);
		y_size_ -= (nTop + nBottom);

		m_oEnvp.MinX += nLeft*m_dblRes;
		m_oEnvp.MaxX -= nRight*m_dblRes;
		m_oEnvp.MinY += nBottom*m_dblRes;
		m_oEnvp.MaxY -= nTop*m_dblRes;
		
		
		return true;
	}
  

	ClassifiedRasterBuffer* GeoRasterBuffer::ClassifyByPredefinedIntervals(int nNumClass, int* panIntervals)
	{
		if (nNumClass == 0 || panIntervals == 0) return 0;
		
		int nArea = x_size_*y_size_;
		unsigned char* panClasses = new unsigned char[nArea];
		int16_t* panValues = (int16_t*)p_pixel_data_;
				
		int c;

		for (int i = 0; i < nArea; i++)
		{
      if (panValues[i]==m_nNODV) panClasses[i] = 0;
      else
      {
        for (c = 1; c <= nNumClass; c++)
          if (panIntervals[c]>panValues[i]) break;
        panClasses[i] = c <= nNumClass ? c : 0;
      }
			
		}

		ClassifiedRasterBuffer* poClassifiedBuffer = new ClassifiedRasterBuffer(nNumClass);
		poClassifiedBuffer->CreateBuffer(1, x_size_, y_size_, panClasses);
		poClassifiedBuffer->CloneGeoRef(this);

		delete[]panClasses;
		return poClassifiedBuffer;
	}

	
	double* GeoRasterBuffer::CalcPixelRanks(int nBand)
	{
		int16_t* panPixels = (int16_t*)p_pixel_data_;
		int nArea = x_size_ * y_size_;
		int nNumValidPixels = 0;
		for (int i = 0; i < nArea; i++)
		{
			if (panPixels[i] != m_nNODV) nNumValidPixels++;
		}
		if (nNumValidPixels == 0) return 0;

		double* padblPixelRanks = new double[nArea];
		for (int i = 0; i < nArea; i++)
			padblPixelRanks[i] = 0.;

		for (int b = 0; b < num_bands_; b++)
		{
			int nMax = -0xfffffff, nMin = 0xfffffff;
			int nBandOffset = b * nArea;
			for (int i = nBandOffset; i < nBandOffset + nArea; i++)
			{
				if (panPixels[i] != m_nNODV)
				{
					nMax = nMax < panPixels[i] ? panPixels[i] : nMax;
					nMin = nMin > panPixels[i] ? panPixels[i] : nMin;
				}
			}

			int nHistSize = nMax - nMin + 1;
			int* panHist = new int[nHistSize];
			for (int i = 0; i < nHistSize; i++)
				panHist[i] = 0;

			for (int i = nBandOffset; i < nBandOffset + nArea; i++)
			{
				if (panPixels[i] != m_nNODV) panHist[panPixels[i] - nMin]++;
			}

			for (int i = 1; i < nHistSize; i++)
				panHist[i] += panHist[i - 1];

			for (int i = 0; i < nArea; i++)
			{
				if (panPixels[nBandOffset + i] != m_nNODV)
				{
					padblPixelRanks[i] += ( ((double)panHist[panPixels[nBandOffset + i]-nMin]) 
											/ nNumValidPixels);
				}
			}
			delete[]panHist;
		}

				
		if (num_bands_ > 1)
		{
			for (int i = 0; i < nArea; i++)
			{
				padblPixelRanks[i] /= num_bands_;
			}
		}
		
		return padblPixelRanks;
	}

	ClassifiedRasterBuffer* GeoRasterBuffer::ClassifyByPredefinedQuantiles(int nNumClass, 
																			double* padblQuantiles)
	{
		double* padblPixelRanks = CalcPixelRanks();
		if (!padblPixelRanks) return 0;

		int nArea = x_size_ * y_size_;
		unsigned char* panClasses = new unsigned char[nArea];

		for (int i = 0; i < nArea; i++)
		{
			for (int n = 0; n <= nNumClass; n++)
			{
				if (padblPixelRanks[i] < padblQuantiles[n] + 1e-20)
				{
					panClasses[i] = n;
					break;
				}
			}
		}

		ClassifiedRasterBuffer* poClassifiedBuffer = new ClassifiedRasterBuffer(nNumClass);
		poClassifiedBuffer->CreateBuffer(1, x_size_, y_size_, panClasses);
		poClassifiedBuffer->CloneGeoRef(this);

		delete[]panClasses;
		delete[]padblPixelRanks;
		return poClassifiedBuffer;
		/*
		int16_t* panPixels = (int16_t*)p_pixel_data_;
		int nArea = x_size_ * y_size_;
		int nNumPixels = 0;
		int nMax = -0xfffffff, nMin = 0xfffffff;
		for (int i = 0; i < nArea; i++)
		{
			if (panPixels[i] != m_nNODV)
			{
				nNumPixels++;
				nMax = nMax < panPixels[i] ? panPixels[i] : nMax;
				nMin = nMin > panPixels[i] ? panPixels[i] : nMin;
			}
		}

		if (nNumPixels == 0)
		{
			panIntervals = 0;
			return 0;
		}

		unsigned char* panClasses = new unsigned char[nArea];

		for (int b = 0; b < this->num_bands_; b++)
		{
			int nHistSize = nMax - nMin + 1;
			int* panHist = new int[nHistSize];
			for (int i = 0; i < nHistSize; i++)
				panHist[i] = 0;

			for (int i = 0; i < nArea; i++)
			{
				if (panPixels[i] != m_nNODV) panHist[panPixels[i] - nMin]++;
			}

			panIntervals = new int[nNumClass + 1];
			panIntervals[0] = nMin;
			panIntervals[nNumClass] = nMax + 1;

			int nClassSum = 0;
			int n = 1;

			for (int i = 0; i < nHistSize; i++)
			{
				if ((((double)(nClassSum + panHist[i])) / ((double)nNumPixels)) + padblQuantiles[n - 1] > padblQuantiles[n])
				{
					panIntervals[n] = i + nMin;
					nClassSum = 0;
					if (n == nNumClass - 1) break;
					else n++;
				}
				nClassSum += panHist[i];
			}


			if (n < nNumClass - 1)
			{
				for (int i = n; i++; i < nNumClass)
					panIntervals[i] = nMax;
			}
		}
		delete[]panHist;
		return ClassifyByPredefinedIntervals(nNumClass, panIntervals);
		*/
	}

	//ToDo
	ClassifiedRasterBuffer* GeoRasterBuffer::ClassifyEqualArea(int nNumClass, int* &panIntervals)
	{
		int16_t* panPixels = (int16_t*)p_pixel_data_;
		int nArea = x_size_*y_size_;
		int nNumPixels = 0;
		int nMax = -0xfffffff, nMin = 0xfffffff;
		for (int i = 0; i < nArea; i++)
		{
			if (panPixels[i]!=m_nNODV)
			{
				nNumPixels++;
				nMax = nMax < panPixels[i] ? panPixels[i] : nMax;
				nMin = nMin > panPixels[i] ? panPixels[i] : nMin;
			}
		}

		if (nNumPixels == 0)
		{
			panIntervals = 0;
			return 0;
		}

		int nHistSize = nMax - nMin + 1;
		int* panHist = new int[nHistSize];
		for (int i = 0; i < nHistSize; i++)
			panHist[i] = 0;

		for (int i = 0; i < nArea; i++)
		{
			if (panPixels[i]!=m_nNODV) panHist[panPixels[i] - nMin]++;
		}

		panIntervals = new int[nNumClass + 1];
		panIntervals[0] = nMin;
		panIntervals[nNumClass] = nMax + 1;

		int nFreqSum = 0;
		int nClassPixelCount = nNumPixels / nNumClass;
		int n = 1;

		
		for (int i = 0; i < nHistSize; i++)
		{
			if ((nFreqSum / nClassPixelCount) < ((nFreqSum + panHist[i]) / nClassPixelCount))
			{
				panIntervals[n] = i + nMin;
				if (n == nNumClass - 1) break;
				else n++;
			}
			nFreqSum += panHist[i];
		}
		
		if (n < nNumClass-1)
		{
			for (int i = n; i++; i < nNumClass)
				panIntervals[i] = nMax;
		}

		delete[]panHist;
		return ClassifyByPredefinedIntervals(nNumClass, panIntervals);
	}

	//TODO	- case dblSTDCoeff = 0 and write values to panRanges
	ClassifiedRasterBuffer* GeoRasterBuffer::ClassifyEqualIntervals(int nNumClass,
		int* &panIntervals, double dblSTDCoeff)
	{
		int16_t* panPixels = (int16_t*)p_pixel_data_;

		double dblMean = 0;
		double dblSTD = 0;
		int nNumPixels = 0;
		int nMax = -0xfffffff, nMin = 0xfffffff;

		int nArea = x_size_*y_size_;
		for (int i = 0; i < nArea; i++)
		{
			if (panPixels[i]!=m_nNODV)
			{
				nNumPixels++;
				if (dblSTDCoeff > 0)
				{
					dblMean += panPixels[i];
					dblSTD += panPixels[i] * panPixels[i];
				}
				nMax = nMax < panPixels[i] ? panPixels[i] : nMax;
				nMin = nMin > panPixels[i] ? panPixels[i] : nMin;				
			}
		}
		
		if (nNumPixels == 0)
		{
			panIntervals = 0;
			return 0;
		}

		panIntervals = new int[nNumClass + 1];
		panIntervals[0] = nMin;
		panIntervals[nNumClass] = nMax + 1;

		if (dblSTDCoeff > 0)
		{
			dblMean /= nNumPixels;
			dblSTD = sqrt((dblSTD / nNumPixels) - dblMean*dblMean);
			double dblMin = dblMean - dblSTD * dblSTDCoeff;
			double dblMax = dblMean + dblSTD * dblSTDCoeff;
			double dblDiff = (dblMax - dblMin) / nNumClass;
			for (int i = 1; i < nNumClass; i++)
				panIntervals[i] = (int)(dblMin + i*dblDiff + 0.5);
		}

		return ClassifyByPredefinedIntervals(nNumClass, panIntervals);
	}

	
	
	GeoRasterBuffer* GeoRasterBuffer::BurnVectorMask(OGRGeometry* poClipGeom, double dblPixelOffset)
	{
		GeoRasterBuffer* poRasterizedVector = new GeoRasterBuffer();
		poRasterizedVector->CreateBuffer(1, x_size_, y_size_);
		poRasterizedVector->InitByValue(1);
		poRasterizedVector->CloneGeoRef(this);

		OGRGeometry* poBufferedClipGeom = dblPixelOffset ?
			poClipGeom->Buffer(dblPixelOffset*this->m_dblRes) : poClipGeom;

		if (!poBufferedClipGeom)
		{
			delete(poRasterizedVector);
			return 0;
		}
		
		bool bClipResult = poRasterizedVector->Clip(poBufferedClipGeom);

		if (dblPixelOffset) delete(poBufferedClipGeom);

		if (! bClipResult)
		{
			delete (poRasterizedVector);
			return 0;
		}
		else return poRasterizedVector;
	}


 	bool GeoRasterBuffer::Clip(OGRGeometry* poClipGeom, double dblPixelOffset)
	{
				
		double padblGeoTransform[6];
		if (!CalcGeoTransform(padblGeoTransform)) return 0;
				
		GDALDataset* poBufferDS = CreateInMemGDALDataset();
		GDALDataset* poClippedDS = CreateInMemGDALDataset(false);

		GDALWarpOptions *poWarpOptions = GDALCreateWarpOptions();
		poWarpOptions->papszWarpOptions = 0;
		poWarpOptions->hSrcDS = poBufferDS;
		poWarpOptions->hDstDS = poClippedDS;
		poWarpOptions->dfWarpMemoryLimit = 250000000;
		double			dblError = 0.125;
		
		char* pachSRSWKT;
		m_poSRS->exportToWkt(&pachSRSWKT);
		
		OGRGeometry* poBufferClipGeom = dblPixelOffset == 0 ?
			poClipGeom : poClipGeom->Buffer(this->m_dblRes*dblPixelOffset);
		poWarpOptions->hCutline = VectorOperations::ConvertFromSRSToPixelLine(poBufferClipGeom, padblGeoTransform);

		poWarpOptions->pTransformerArg = GDALCreateApproxTransformer(GDALGenImgProjTransform,
			GDALCreateGenImgProjTransformer(poBufferDS, pachSRSWKT, poClippedDS, pachSRSWKT, false, 0.0, 1),
			dblError);
		poWarpOptions->pfnTransformer = GDALApproxTransform;
		//poWarpOptions->eResampleAlg = GRA_NearestNeighbour;
		poWarpOptions->eResampleAlg = GRA_Cubic;


		
		poWarpOptions->padfSrcNoDataReal = new double[num_bands_];
		poWarpOptions->padfSrcNoDataImag = new double[num_bands_];
		for (int i = 0; i < num_bands_; i++)
		{
			poWarpOptions->padfSrcNoDataReal[i] = m_nNODV;
			poWarpOptions->padfSrcNoDataImag[i] = 0;
		}

		// Initialize and execute the warp operation. 
		GDALWarpOperation gdal_warp_operation;
		gdal_warp_operation.Initialize(poWarpOptions);

		bool  warp_error = (CE_None != gdal_warp_operation.ChunkAndWarpImage(0, 0, x_size_, y_size_));

		GDALDestroyApproxTransformer(poWarpOptions->pTransformerArg);

		delete((OGRGeometry*)poWarpOptions->hCutline);
		delete(poWarpOptions->padfSrcNoDataReal);
		delete(poWarpOptions->padfSrcNoDataImag);
		poWarpOptions->padfSrcNoDataReal = 0;
		poWarpOptions->padfSrcNoDataImag = 0;

		poWarpOptions->hCutline = 0;
		GDALDestroyWarpOptions(poWarpOptions);
		OGRFree(pachSRSWKT);

		poClippedDS->RasterIO(GF_Read, 0, 0, x_size_, y_size_, get_pixel_data_ref(),
			x_size_, y_size_, data_type_, num_bands_, 0, 0, 0, 0);
		GDALClose(poBufferDS);
		GDALClose(poClippedDS);
		
		//ToDoo VSIUnlink: poBufferDS, poClippedDS - need testing

	
		if (dblPixelOffset != 0) delete(poBufferClipGeom);

		return true;

	}

	
	bool GeoRasterBuffer::SaveGeoRefFile(string strRasterFile)
	{
		SaveBufferToFile(strRasterFile);
		if ((strRasterFile.find(".tif") != string::npos) ||
			(strRasterFile.find(".TIF") != string::npos))
		{
			GDALDataset* poDS = (GDALDataset*)GDALOpen(strRasterFile.c_str(), GA_Update);
			
			char* pachWKT;
			GetSRSRef()->exportToWkt(&pachWKT);
			poDS->SetProjection(pachWKT);
			
			double dblGeotransform[6];
			dblGeotransform[0] = m_oEnvp.MinX;
			dblGeotransform[1] = m_dblRes;
			dblGeotransform[2] = 0;
			dblGeotransform[3] = m_oEnvp.MaxY;
			dblGeotransform[4] = 0;
			dblGeotransform[5] = -m_dblRes;
			poDS->SetGeoTransform(dblGeotransform);
			
			OGRFree(pachWKT);
			GDALClose(poDS);
			//char *p_dst_wkt = 0;
			//p_tile_mset_->GetTilingSRSRef()->exportToWkt(&p_dst_wkt);
			//p_vrt_ds->SetProjection(p_dst_wkt);
			//OGRFree(p_src_wkt);

		}
		else GMXFileSys::WriteWLDFile(strRasterFile, m_oEnvp.MinX, m_oEnvp.MaxY, m_dblRes);
		
		return true;
	}


	map<int, OGRMultiPolygon*> GeoRasterBuffer::Polygonize(GeoRasterBuffer* poGeoBuffer)
	{
		map<int, OGRMultiPolygon*> mapOutput;

		GDALDataset* poInMemDS = poGeoBuffer->CreateInMemGDALDataset();
				
		string	strInMemName = ("/vsimem/shpinmem_" + GMXString::ConvertIntToString(rand()));
		GDALDataset* poInMemSHP = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile")->
			Create(strInMemName.c_str(), 0, 0, 0, GDT_Unknown, NULL);
		OGRLayer* poLayer = poInMemSHP->CreateLayer(strInMemName.c_str(), m_poSRS, wkbPolygon, NULL);
		OGRFieldDefn oFieldDefn("DN", OFTInteger);
		poLayer->CreateField(&oFieldDefn);


		if ((OGRERR_NONE == GDALPolygonize(poInMemDS->GetRasterBand(1), poInMemDS->GetRasterBand(1),
			poLayer, 0, 0, 0, 0)) && (poLayer->GetFeatureCount()>0))
		{

			while (OGRFeature* poFeature = poLayer->GetNextFeature())
			{
				int dn = poFeature->GetFieldAsInteger("DN");
				if (mapOutput.find(dn) == mapOutput.end())
					mapOutput[dn] =
					(OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);

				mapOutput[dn]->addGeometry(poFeature->GetGeometryRef());

				OGRFeature::DestroyFeature(poFeature);
			}
		}


		GDALClose(poInMemDS);
		poInMemSHP->DeleteLayer(0);
		GDALClose(poInMemSHP);
		VSIUnlink(strInMemName.c_str());
		return mapOutput;

	}

	map<int, OGRMultiPolygon*> GeoRasterBuffer::Polygonize()
	{
		return Polygonize(this);
	}

	map<int,OGRMultiPolygon*> GeoRasterBuffer::Polygonize(int nOffset, int nStep, int nUpperBound)
	{
			
		int nNumClass = 1 + (nUpperBound - nOffset) / nStep;
		int *panIntervals = new int[nNumClass + 1];
		for (int i = 0; i <= nNumClass; i++)
			panIntervals[i] = nOffset + i*nStep;
		

		ClassifiedRasterBuffer* poRankedValuesBuffer = ClassifyByPredefinedIntervals(nNumClass,panIntervals);
		
		map<int, OGRMultiPolygon*> mapOutputTemp = Polygonize(poRankedValuesBuffer);
		delete(poRankedValuesBuffer);

		map<int, OGRMultiPolygon*> mapOutputFinal;

		for (const auto &oPair : mapOutputTemp)
			mapOutputFinal[panIntervals[oPair.first]] = oPair.second;
		
		return mapOutputFinal;
	}

	/*
	OGRMultiPolygon* GeoRasterBuffer::PolygonizeRange(int nMinValue, int nMaxValue)
	{
	
		OGRMultiPolygon** papoMPolygons = Polygonize(nMinValue, nMaxValue - nMinValue + 1, nMaxValue);
		OGRMultiPolygon* poM

	
		GeoRasterBuffer* poMaskBuffer = 0;
		switch (data_type_)
		{
			case GDT_Byte:
			{
				unsigned char t = 1;
				poMaskBuffer = CreateMaskBufferByRange(t, nMinValue, nMaxValue);
				break;
			}
			case GDT_UInt16:
			{
				uint16_t t = 257;
				poMaskBuffer = CreateMaskBufferByRange(t, nMinValue, nMaxValue);
				break;
			}
			case GDT_Int16:
			{
				int16_t t = 257;
				poMaskBuffer = CreateMaskBufferByRange(t, nMinValue, nMaxValue);
				break;
			}
			default:
				return 0;
		}
		GDALDataset* poInMemDS = poMaskBuffer->CreateInMemGDALDataset(true);

		string	strInMemName = ("/vsimem/shpinmem_" + GMXString::ConvertIntToString(rand()));
		GDALDataset* poInMemSHP = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile")->
			Create(strInMemName.c_str(), 0, 0, 0, GDT_Unknown, NULL);
		OGRLayer* poLayer = poInMemSHP->CreateLayer(strInMemName.c_str(), m_poSRS, wkbPolygon, NULL);
		OGRFieldDefn oFieldDefn("DN", OFTInteger);
		poLayer->CreateField(&oFieldDefn);
		
		OGRMultiPolygon* poOutputMPoly = 0;
		
		if ((OGRERR_NONE == GDALPolygonize(poInMemDS->GetRasterBand(1), poInMemDS->GetRasterBand(1),
			poLayer, 0, 0, 0, 0)) && (poLayer->GetFeatureCount()>0))
		{
			poOutputMPoly = new OGRMultiPolygon();
			while (OGRFeature* poFeature = poLayer->GetNextFeature())
			{
				poOutputMPoly->addGeometry(poFeature->GetGeometryRef());
				OGRFeature::DestroyFeature(poFeature);
			}
		}
		
		GDALClose(poInMemDS);
		poInMemSHP->DeleteLayer(0);
		GDALClose(poInMemSHP);
		VSIUnlink(strInMemName.c_str());

		return poOutputMPoly;
		
	}
	*/
	
	/*
	bool GeoRasterBuffer::TraceEdge(int nMinValue,
		int nMaxValue,
		int nX0,
		int nY0,
		int nX1,
		int nY1,
		list<pair<int, int>> &listEdgeLine)
	{
	
		uint16_t* panValues = (uint16_t*)this->p_pixel_data_;

		int panArrows[4][2] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
		
		listEdgeLine.push_back(pair<int, int>(nX0, nY0));
		listEdgeLine.push_back(pair<int, int>(nX1, nY1));
	
			
		while ( ((*listEdgeLine.begin()).first != listEdgeLine.back().first) || 
			((*listEdgeLine.begin()).second != listEdgeLine.back().second))
		{
			for (int i = 0; i < 4; i++)
			{
				if ((nX1 + panArrows[i][0] == nX0) && (nY1 + panArrows[i][1] == nY0)) continue;
			
				int nW1 = nX1 - 1 + (panArrows[i][0]>0);
				int nW2 = nX1 - (panArrows[i][0]<0);
				
				int nH1 = nY1 - 1 + (panArrows[i][1]>0);
				int nH2 = nY1 - (panArrows[i][1]<0);

				bool bPoint1InRange, bPoint2InRange;

				if ((nW1 < 0) || (nW1 == this->x_size_) || (nH1 < 0) || (nH1 == this->y_size_)) bPoint1InRange = false;
				else
				{
					bPoint1InRange = ((panValues[nH1*x_size_ + nW1]<nMinValue) ||
						(panValues[nH1*x_size_ + nW1]>nMaxValue)) ? false : true;
				}

				if ((nW2 < 0) || (nW2 == this->x_size_) || (nH2 < 0) || (nH2 == this->y_size_)) bPoint2InRange = false;
				else
				{
					bPoint2InRange = ((panValues[nH2*x_size_ + nW2]<nMinValue) || 
									(panValues[nH2*x_size_ + nW2]>nMaxValue)) ? false : true;
				}

				if (bPoint1InRange != bPoint2InRange)
				{
					listEdgeLine.push_back(pair<int, int>(nX1 + panArrows[i][0], nY1 + panArrows[i][1]));
					nX0 = nX1;
					nY0 = nY1;
					nX1 = nX1 + panArrows[i][0];
					nY1 = nY1 + panArrows[i][1];
					break;
				}
			}
			//choose from one or two edges what is next
		}



		return true;
	}
	*/

	/*
	bool GeoRasterBuffer::CalculateContour(int nMinValue, int nMaxValue, OGRGeometry* &poSRSGeometry)
	{
		//iterate all pixels
		//if pixel is on edge and there are untraced sides on edge
		//pick a side on edge and start tracing
		//repeat for other sides
		//create an array of list<pair<int,int>>
	}
	*/

	/*
	bool Classifier::ConvertPixelsToPolygons(int nZ, string strNDVI, string strOutputVector)
	{
		Classifier oClassifier;
		if (!oClassifier.Init(listTileContainers,
			oOptionParser.GetOptionValue("-v"),
			oOptionParser.GetOptionValue("-z") == "" ? 0 : atoi(oOptionParser.GetOptionValue("-z").c_str())))
		{
			cout << "ERROR: reading input tile containers" << endl;
			return 5;
		}

		//debug
		oClassifier.GetGeoRasterBufferRef()->Polygonize("F:\\mpotanin\\data6\\poly_174_z12.shp");

	}
	*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ClassifiedRasterBuffer::PolygonizePixels(string strOutputVectorFile,
                                              OGRGeometry* poClipMask,
                                              bool bSaveTo4326)
{
  ClassifiedRasterBuffer* poClippedBuffer = new ClassifiedRasterBuffer(this->m_nClasses);
  poClippedBuffer->CreateBuffer(this);
  poClippedBuffer->CloneGeoRef(this);

  poClippedBuffer->Clip(poClipMask);
  poClippedBuffer->ReplaceByInterpolatedValues(poClipMask, 1, 1);
  poClippedBuffer->AdjustExtentToClippedArea();
  
  //debug
  //poClippedBuffer->SaveGeoRefFile(strOutputVectorFile+"_2.tif");
  //end-debug

  GDALDriver *poDriver = GMXFileSys::GetExtension(strOutputVectorFile) == "shp" ?
    GetGDALDriverManager()->GetDriverByName("ESRI Shapefile") :
    GetGDALDriverManager()->GetDriverByName("GeoJSON");
  GDALDataset *poDS = poDriver->Create(strOutputVectorFile.c_str(), 0, 0, 0, GDT_Unknown, NULL);
  if (poDS == NULL)
  {
    printf("Creation of output file failed.\n");
    exit(1);
  }

  OGRSpatialReference oSRS4326;
  oSRS4326.SetWellKnownGeogCS("WGS84");

  OGRLayer *poLayer = poDS->CreateLayer("pixels", bSaveTo4326 ? &oSRS4326 : this->m_poSRS, wkbPolygon, NULL);
  OGRFieldDefn oFieldDN("DN", OFTInteger);
  poLayer->CreateField(&oFieldDN);
  OGRFieldDefn oFieldFRAC("FRACTION", OFTReal);
  poLayer->CreateField(&oFieldFRAC);

    
  unsigned char* p_class_values = (unsigned char*)p_pixel_data_;

  int val = 0;
  int n = x_size_*y_size_;
  for (int i = 0; i < x_size_; i++)
  {
    for (int j = 0; j < y_size_; j++)
    {
      val = p_class_values[x_size_*j + i];
      if (val == 0) continue;
      OGRLinearRing oPixelRing;
      oPixelRing.addPoint(m_oEnvp.MinX + i*m_dblRes, m_oEnvp.MaxY - j*m_dblRes);
      oPixelRing.addPoint(m_oEnvp.MinX + (i + 1)*m_dblRes, m_oEnvp.MaxY - j*m_dblRes);
      oPixelRing.addPoint(m_oEnvp.MinX + (i + 1)*m_dblRes, m_oEnvp.MaxY - (j + 1)*m_dblRes);
      oPixelRing.addPoint(m_oEnvp.MinX + i*m_dblRes, m_oEnvp.MaxY - (j + 1)*m_dblRes);
      oPixelRing.closeRings();
      OGRPolygon	*poPixelPoly = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
      poPixelPoly->addRing(&oPixelRing);
      
      if (!poClipMask->Intersects(poPixelPoly))
      {
        delete(poPixelPoly);
        continue;
      }
      else
      {
        OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        poFeature->SetField("DN", val);
       
        if (poClipMask->Contains(poPixelPoly))
          poFeature->SetField("FRACTION", 1.0);
        else
        {
          OGRGeometry* poPixelCut = poClipMask->Intersection(poPixelPoly);
          if (!poPixelCut)
            poFeature->SetField("FRACTION", 1.0);
          else
          {
            double dblCutArea = 0.0;
            switch (poPixelCut->getGeometryType())
            {
              case wkbPolygon:
                dblCutArea = ((OGRPolygon*)poPixelCut)->get_Area();
                break;
              case wkbMultiPolygon:
                dblCutArea = ((OGRMultiPolygon*)poPixelCut)->get_Area();
                break;
              default:
                dblCutArea = poPixelPoly->get_Area();
                break;
            }
            delete(poPixelCut);
            poFeature->SetField("FRACTION", dblCutArea / poPixelPoly->get_Area());
         }
        }

        if (bSaveTo4326)
        {
          poPixelPoly->assignSpatialReference(this->m_poSRS);
          poPixelPoly->transformTo(&oSRS4326);
          poPixelPoly->assignSpatialReference(0);
        }

        poFeature->SetGeometryDirectly(poPixelPoly);
        poLayer->CreateFeature(poFeature);
        OGRFeature::DestroyFeature(poFeature);
      }
    }
  }

  GDALClose(poDS);
  delete(poClippedBuffer);
  return true;
}





	bool ClassifiedRasterBuffer::ReplaceByInterpolatedValues(OGRGeometry* poVector,
		double dblPixelInward, double dblPixelOutward)
	{
		GeoRasterBuffer* poMaskOutwardOffset = this->BurnVectorMask(poVector, dblPixelOutward);
		GeoRasterBuffer* poMaskInwardOffset = this->BurnVectorMask(poVector, -dblPixelInward);

		unsigned char* panOutwardBufferPixels = 
			(unsigned char*)poMaskOutwardOffset->get_pixel_data_ref();
		
		unsigned char* panInwardBufferPixels =
			(unsigned char*)poMaskInwardOffset->get_pixel_data_ref();

		int nArea = x_size_*y_size_;

		for (int i = 0; i < nArea; i++)
		{
			panOutwardBufferPixels[i] = panOutwardBufferPixels[i] ^ panInwardBufferPixels[i];
		}
	
    //debug
    //poMaskOutwardOffset->SaveGeoRefFile("F:\\Work\\Projects\\agri-zoning\\testdata\\task5_2\\border.tif");
		//end-debug

    ApplyMajorityFilter(int(2 * (dblPixelOutward + dblPixelInward) + 3.5),true,
			                  panOutwardBufferPixels);

		delete(poMaskInwardOffset);
		delete(poMaskOutwardOffset);


		return true;
	}

	
	bool ClassifiedRasterBuffer::ApplyMajorityFilter(int nWinSize, 
                                                    bool bOnlyNoDataPixels,
                                                    unsigned char* pabInterpolationMask )
	{
		int *panFreq = new int[m_nClasses + 1];
		int nWinSize_2 = nWinSize / 2;
		int nMaxFreqClass;

		int nArea = x_size_*y_size_;
		unsigned char* pabClasses = (unsigned char*)p_pixel_data_;
		unsigned char* pabClassesFiltered = new unsigned char[nArea];
		for (int i = 0; i<nArea; i++)
			pabClassesFiltered[i] = pabClasses[i];

		
		for (int y = 0; y < y_size_; y++)
		{
			for (int x = 0; x < x_size_; x++)
			{
        if ( bOnlyNoDataPixels && (pabClasses[y*x_size_ + x] != 0)) continue;
          
				if (!pabInterpolationMask)
				{
					if (pabClasses[y*x_size_ + x] == 0) continue;
				}
				else
				{
					if (pabInterpolationMask[y*x_size_ + x]==0) continue;
				}
								
				for (int q = 0; q <= m_nClasses; q++)
					panFreq[q] = 0;
				
				for (int i = ((y>=nWinSize_2) ? -nWinSize_2 : -y); 
						 i <= ((y<y_size_-nWinSize_2) ? nWinSize_2 : y_size_-y-1);
						 i++)
				{
					for (int j = ((x >= nWinSize_2) ? -nWinSize_2 : -x);
						j <= ((x < x_size_ - nWinSize_2) ? nWinSize_2 : x_size_ - x - 1);
						j++)
					{
						if (!pabInterpolationMask)
							panFreq[pabClasses[x_size_*(y + i) + x + j]]++;
						else if (!pabInterpolationMask[x_size_*(y + i) + x + j])
							panFreq[pabClasses[x_size_*(y + i) + x + j]]++;
					}
					
				}

				nMaxFreqClass = pabClasses[x_size_*y + x];
				panFreq[0] = 0;
				for (int q = 1; q <= m_nClasses; q++)
				{
        	nMaxFreqClass = panFreq[q] > panFreq[nMaxFreqClass] ? q : nMaxFreqClass;
        }
				pabClassesFiltered[x_size_*y + x] = nMaxFreqClass;
			}
		}

		this->p_pixel_data_ = pabClassesFiltered;
		
		delete[]pabClasses;
		delete[]panFreq;

		return true;
	}

	OGREnvelope GeoRasterBufferCollection::CalculateBundleEnvelope()
	{
		OGREnvelope oBundleEnvp = (*m_mapBuffers.begin()).second->GetEnvelope();
		for (auto iter : m_mapBuffers)
		{
			oBundleEnvp=gmx::VectorOperations::MergeEnvelopes(oBundleEnvp,
				iter.second->GetEnvelope());
		}
		return oBundleEnvp;
	}

	GeoRasterBufferCollection::~GeoRasterBufferCollection()
	{
		for (auto iter : m_mapBuffers)
			delete(iter.second);
	}

	GeoRasterBufferCollection* InitBundleFromFileInput(string strNDVIFilesTable,
		string strVectorFile, string strBasePath, int nZoom = 13)
	{

		return 0;
	}



	bool GeoRasterBufferCollection::Init(map<string, list<string>> mapNDVITiles,
		map<string, OGRGeometry*> mapClipGeom, int nZoom)
	{

		for (auto iter : mapNDVITiles)
		{
			list<string> &listNDVI = iter.second;
			string strName = iter.first;
			if (!listNDVI.size() || !mapClipGeom[strName])
			{
				cout << "ERROR: null value input parameters for object: " << strName << endl;
				continue;
			}
			else
			{
				m_mapBuffers[strName] =
					GeoRasterBuffer::InitFromNDVITiles(listNDVI,mapClipGeom[strName],nZoom);
				if (!m_mapBuffers[strName])
				{
					cout << "ERROR: can't init properly object: " << strName << endl;
					//ToDo
					//Clear();					
					return false;
				}
			}
		}
		return true;
	}

	bool GeoRasterBufferCollection::SaveToFile(string strFile, GDALColorTable *poColTab)
	{
		GeoRasterBuffer* poSampleBuffer = (*m_mapBuffers.begin()).second;

		OGREnvelope oBundleEnvp = CalculateBundleEnvelope();
		int nBundleW = (((oBundleEnvp.MaxX - oBundleEnvp.MinX) /
			poSampleBuffer->GetPixelSize()) + 0.5);
		int nBundleH = (((oBundleEnvp.MaxY - oBundleEnvp.MinY) /
			poSampleBuffer->GetPixelSize()) + 0.5);

		GeoRasterBuffer oBundleBuffer;
		oBundleBuffer.CreateBuffer(poSampleBuffer->get_num_bands(), nBundleW, nBundleH, 
			0, poSampleBuffer->get_data_type(),0,poColTab);
		
		for (auto iter : m_mapBuffers)
		{
			GeoRasterBuffer* poIterBuffer = iter.second;
			int nOffsetLeft = ((poIterBuffer->GetEnvelope().MinX - oBundleEnvp.MinX) /
				poSampleBuffer->GetPixelSize()) + 0.5;

			int nOffsetTop = ((oBundleEnvp.MaxY - poIterBuffer->GetEnvelope().MaxY) /
				poSampleBuffer->GetPixelSize()) + 0.5;
			oBundleBuffer.SetPixelDataBlock(nOffsetLeft, nOffsetTop, poIterBuffer->get_x_size(),
				poIterBuffer->get_y_size(), poIterBuffer->get_pixel_data_ref());
		}
			
		oBundleBuffer.SaveBufferToFile(strFile);

		return true;
	}
	/*
	GeoRasterBuffer* Classifier::ClassifyByPredefinedIntervals(int nNumClass, int* panRanges)
	{
		int nWidth = m_poGeoBuffer->get_x_size();
		int nHeight = m_poGeoBuffer->get_y_size();
		uint16_t* panIVI = CalcIVIRaster();

		int nBands = m_poGeoBuffer->get_num_bands();

		unsigned char* panClasses = new unsigned char[nWidth*nHeight];
		int t;
		for (int i = 0; i < nHeight; i++)
		{
			for (int j = 0; j < nWidth; j++)
			{
				t = i*nWidth + j;
				if (panIVI[t] < panRanges[0]) panClasses[t] = 0;
				else if (panIVI[t] >= panRanges[nNumClass]) panClasses[t] = 0;
				else
				{
					for (int c = 0; c < nNumClass; c++)
					{
						if (panIVI[t] < panRanges[c + 1])
						{
							panClasses[t] = c + 1;
							break;
						}
					}
				}
			}
		}
		
		GeoRasterBuffer *poOutputBuffer = new GeoRasterBuffer();
		poOutputBuffer->CreateBuffer(1, nWidth, nHeight, panClasses);

		poOutputBuffer->CloneGeoRef(m_poGeoBuffer);


		delete[]panClasses;
		delete[]panIVI;

		return poOutputBuffer;

	}
	*/
	

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
			if (sMeta.nYear==2013 || sMeta.nYear == nLastYear) continue;
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

	string NDVIProfile::GetFullPathBySceneid(string strSceneid, string strBasePath)
	{
		//\\tinkerbell - smb\ifs\kosmosnimki\Operative\alt_proc\ls8\ndvi\2014\2014 - 03 - 19\LC81720282014078LGN00_ndvi.tiles
		//	\\tinkerbell - smb\ifs\kosmosnimki\Operative\alt_proc\s2\ndvi\2016\2016 - 06 - 19\S2A_L1C_20160619_124_ndvi.tiles
		//	\\tinkerbell - smb\ifs\kosmosnimki\Operative\alt_proc\ls8\ndvi\2017\2017 - 04 - 12\LC81720282017102LGN00_ndvi.tiles
		strBasePath += strSceneid.find("LC8", 0) != string::npos ? "\\ls8\\ndvi\\" :
																	"\\s2\\ndvi\\";

		return "";
	}


	bool NDVIProfile::ParseDateString(string strDate, int* pnYear, int* pnMonth, int* pnDay, int* pnDOY)
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


	bool NDVIProfile::ParseInputData(string strHRMeanFile, 
									string strMODISMaxFile,
									string strAliasPath, 
									string strRealPath)
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
			
			oMeta.strFileName = strLine.substr(19, (strLine.find(';', 19) - 19));
			if (strAliasPath != "")
			{
				if (oMeta.strFileName.find(strAliasPath, 0) == 0)
					oMeta.strFileName = strRealPath + oMeta.strFileName.substr(strAliasPath.size());
			}
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
			double dblMax = atof(strLine.substr(0, strLine.find(";")).replace(1,1,".").c_str());
			strLine = strLine.substr(strLine.find(";") + 1);
			if (!NDVIProfile::ParseDateString(strLine.substr(0, 10), &nYear, &nMonth, &nDay, &nDOY))
			{
				m_listMetadata.clear();
				m_mapMODISMax.clear();
				return false;
			}
			m_mapMODISMax[nYear] = pair<int,double>(nDOY,dblMax);
		}

		return true;
	}

	/*
	bool Classifier::ClassifyISOCLUSMethod(string strParams, 
											string strBasePath, 
											string strISOFile 
											)
	{
		if (!m_poGeoBuffer) return false;
		unsigned char* pabPixels = (unsigned char*)m_poGeoBuffer->get_pixel_data_ref();
		int nBands = m_poGeoBuffer->get_num_bands();
		int nW = m_poGeoBuffer->get_x_size();
		int nH = m_poGeoBuffer->get_y_size();

		string strFileNameBase = GMXFileSys::GetAbsolutePath(strBasePath, "band");
		for (int i = 0; i < nBands; i++)
			GMXFileSys::SaveDataToFile(strFileNameBase + GMXString::ConvertIntToString(i),
										&pabPixels[i*nW*nH], 
										nW*nH);

		int argc = 3;
		char** argv = new char*[3];
		string strArgs[3] = { "zoning.exe", "-i", "" };
		strArgs[2] = strParams;
		for (int i = 0; i < 3; i++)
		{
			argv[i] = new char[strArgs[i].size() + 1];
			for (int j = 0; j < strArgs[i].size(); j++)
				argv[i][j] = strArgs[i].at(j);
			argv[i][strArgs[i].size()] = 0;
		}

		isoclus_v2::Image::RunISOCLUS(argc, argv);
		for (int i = 0; i < argc; i++)
			delete[]argv[i];
		delete[]argv;
		
		void* pabData;
		int nSize;
		GMXFileSys::ReadBinaryFile(strISOFile, pabData, nSize);
		GeoRasterBuffer oISOBuffer;
		oISOBuffer.CreateBuffer(1, nW, nH, pabData);
		oISOBuffer.CloneGeoRef(m_poGeoBuffer);

		GDALColorTable *poColTab = new GDALColorTable(GPI_RGB);
		GDALColorEntry arrColors[6];
		arrColors[0].c1 = 0; arrColors[0].c2 = 0; arrColors[0].c3 = 0;
		arrColors[1].c1 = 255; arrColors[1].c2 = 0; arrColors[1].c3 = 0;
		arrColors[2].c1 = 247; arrColors[2].c2 = 209; arrColors[2].c3 = 59;
		arrColors[3].c1 = 212; arrColors[3].c2 = 255; arrColors[3].c3 = 190;
		arrColors[4].c1 = 76; arrColors[4].c2 = 227; arrColors[4].c3 = 0;
		arrColors[5].c1 = 47; arrColors[5].c2 = 140; arrColors[5].c3 = 30;
		for (int i = 0; i <= 5; i++)
			poColTab->SetColorEntry(i+1, &arrColors[i]);
		oISOBuffer.set_color_table(poColTab);
		oISOBuffer.SaveGeoRefFile(strISOFile+".png");
		delete(poColTab);

		return true;
	}
	(*/

	bool ZoningMap::SaveToVectorFile(string strFileName, OGRSpatialReference* poSRS)
	{
		const char *pszDriverName = GMXFileSys::GetExtension(GMXString::MakeLower(strFileName)) == "shp" ?
			"ESRI Shapefile" : "GeoJSON";
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
		GDALDataset *poDS = poDriver->Create(strFileName.c_str(), 0, 0, 0, GDT_Unknown, NULL);
		if (!poDS) return false;

		OGRLayer *poLayer = poDS->CreateLayer("zoning", poSRS, wkbPolygon, NULL);
		OGRFieldDefn oField1("NAME", OFTString);
		OGRFieldDefn oField2("DN", OFTInteger);

		poLayer->CreateField(&oField1);
		poLayer->CreateField(&oField2);

		SaveToLayer(poLayer, "DN", "", "");
		
		GDALClose(poDS);

		return true;
	}


	bool ZoningMap::SaveToLayer(OGRLayer* poLayer, string strDNField,
		string strZoningMapNameField, string strZoningMapName)
	{
		if (!poLayer) return false;

		for (auto iter : m_mapZones)
		{
			OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
			poFeature->SetGeometry(iter.second);
			poFeature->SetField(strDNField == "" ? "DN" : strDNField.c_str(), iter.first);
			if (strZoningMapNameField!="")
				poFeature->SetField(strZoningMapNameField.c_str(), strZoningMapName.c_str());
			if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
			{
				OGRFeature::DestroyFeature(poFeature);
				return false;
			}
			else OGRFeature::DestroyFeature(poFeature);
		}

		return true;
	}

	
	bool ZoningMap::TransformToSRS(OGRSpatialReference* poFromSRS, OGRSpatialReference* poToSRS)
	{
		for (auto iter : m_mapZones)
		{
			iter.second->assignSpatialReference(poFromSRS);
			iter.second->transformTo(poToSRS);
			iter.second->assignSpatialReference(0);
		}

		return true;
	}
	
	bool ZoningMap::InitDirectly(map<int, OGRMultiPolygon*> mapZones)
	{
		Clear();
		m_mapZones = mapZones;
		return true;
	}



	bool ZoningMap::InitMakingCopy(map<int, OGRMultiPolygon*> mapZones)
	{
		Clear();
		m_mapZones = mapZones;
		for (auto iter : m_mapZones)
			iter.second = (OGRMultiPolygon*)iter.second->clone();
		return true;
	}
	
	
	bool ZoningMap::Clear()
	{
		for (auto iter : m_mapZones)
			delete(iter.second);
		m_mapZones.clear();
		
		return true;
	}

	map<int, OGRPolygon*> ZoningMap::FindAllBorderingPolygons(OGRPolygon* poPolygon, double dblThreshold)
	{
		map<int, OGRPolygon*> mapBorderingPolygons;

		OGRPolygon* poOneOfPolygons;
		for (auto iter : m_mapZones)
		{
			for (int i = 0; i < iter.second->getNumGeometries(); i++)
			{
				poOneOfPolygons = (OGRPolygon*)iter.second->getGeometryRef(i);
				if (poOneOfPolygons != poPolygon)
				{
					if (poPolygon->Touches(poOneOfPolygons) && 
							(poOneOfPolygons->get_Area()>=dblThreshold))
						mapBorderingPolygons[iter.first] = poOneOfPolygons;
				}
			}
		}

		return mapBorderingPolygons;
	}

	bool ZoningMap::Clip(OGRGeometry* poClipGeometry)
	{
		for (auto iter : m_mapZones)
		{
			OGRGeometry* poIntersection = iter.second->Intersection(poClipGeometry);
			
			if (!poIntersection) return false;
			else
			{
				delete(iter.second);
				m_mapZones[iter.first] = (OGRMultiPolygon*)poIntersection;
			}
		}
		return true;
	}


	bool ZoningMap::FilterByArea(double dblThreshold)
	{
		bool bPolygonsCountDecreased = true;
		while (bPolygonsCountDecreased)
		{
			bPolygonsCountDecreased = false;
			for (auto iter : m_mapZones)
			{
				OGRPolygon* poTesteePolygon;
				for (int i = 0; i < iter.second->getNumGeometries(); i++)
				{
					poTesteePolygon = (OGRPolygon*)iter.second->getGeometryRef(i);
					if (poTesteePolygon->get_Area() < dblThreshold)
					{
						map<int, OGRPolygon*> mapBorderingPolygons = 
							FindAllBorderingPolygons(poTesteePolygon, dblThreshold);

						if (mapBorderingPolygons.size() == 0) continue;

						OGRPolygon* poForUnionPolygon = 0;
						int nNewIndex = 0xFFFFFF;
						double dblLargestArea = 0;

						for (auto iter2 : mapBorderingPolygons)
						{
							if (abs(iter2.first - iter.first) < abs(iter.first - nNewIndex))
							{
								nNewIndex = iter2.first;
								dblLargestArea = iter2.second->get_Area();
								poForUnionPolygon = iter2.second;
							}
							else if (abs(iter2.first - iter.first) == abs(iter.first - nNewIndex))
							{
								if (iter.second->get_Area() > dblLargestArea)
								{
									nNewIndex = iter2.first;
									dblLargestArea = iter2.second->get_Area();
									poForUnionPolygon = iter2.second;
								}
							}
						}
																	
						OGRPolygon* poUnionPolygon = 
							(OGRPolygon*)poForUnionPolygon->Union(poTesteePolygon);
						iter.second->removeGeometry(i);

            gmx::VectorOperations::RemovePolygonFromMultiPolygon(m_mapZones[nNewIndex],
							                                                   poForUnionPolygon);
            
						m_mapZones[nNewIndex]->addGeometryDirectly(poUnionPolygon);
						bPolygonsCountDecreased = true;
						i--;
					}
				}
			}
		}

		RemoveEmptyZones();

		return true;
	}


	bool ZoningMap::RemoveEmptyZones()
	{
		if (m_mapZones.size() == 0) return true; 
		int* panIndexesToDelete = new int[m_mapZones.size()];
		int nCountToDelete = 0;
		for (auto iter : m_mapZones)
		{
			if (iter.second->getNumGeometries() == 0)
			{
				delete(iter.second);
				panIndexesToDelete[nCountToDelete] = iter.first;
				nCountToDelete++;
			}
		}

		for (int i = 0; i < nCountToDelete; i++)
			m_mapZones.erase(panIndexesToDelete[i]);


		return true;
	}



	bool ZoningMapCollection::SaveToVectorFile(string strFileName)
	{
		const char *pszDriverName = GMXFileSys::GetExtension(GMXString::MakeLower(strFileName)) == "shp" ?
			"ESRI Shapefile" : "GeoJSON";
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
		GDALDataset *poDS = poDriver->Create(strFileName.c_str(), 0, 0, 0, GDT_Unknown, NULL);
		if (!poDS) return false;

		OGRLayer *poLayer = poDS->CreateLayer("zoning", m_poSRS, wkbPolygon, NULL);
		OGRFieldDefn oField1("NAME", OFTString);
		OGRFieldDefn oField2("DN", OFTInteger);

		poLayer->CreateField(&oField1);
		poLayer->CreateField(&oField2);

		for (auto iter : m_mapZoningCollection)
		{
			iter.second->SaveToLayer(poLayer, "DN", "NAME", iter.first);
		}

		GDALClose(poDS);
		
		return true;
	}

	map<string, pair<OGRGeometry*, list<string>>>*
		ZoningMapCollectionProcessor::ParseConsoleInput(string strTextFile, string strVectorFile)
	{
		map<string, pair<OGRGeometry*, list<string>>>* pomapParsedInput =
			new map<string, pair<OGRGeometry*, list<string>>>();


		return pomapParsedInput;
	}

}