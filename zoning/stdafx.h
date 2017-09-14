// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#ifdef WIN32
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#define _CRT_SECURE_NO_WARNINGS
#include <tchar.h>
#endif


using namespace std;

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "vrtdataset.h"
#include "gdalwarper.h"
#include "gd.h"
#include "sqlite3.h"
#include "cpl_string.h"
#include "commonutils.h"
#include "gdal_utils_priv.h"


#include "filesystemfuncs.h"
#include "tilingfuncs.h"
#include "tilingparameters.h"
#include "consoleutils.h"





// TODO: reference additional headers your program requires here
