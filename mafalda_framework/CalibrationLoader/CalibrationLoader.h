 /*
 * 	Copyright 2011 John Idarraga, Mathieu Benoit
 *
 * 	This file is part of MAFalda.

    MAFalda is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    MAFalda is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAFalda.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __CalibrationLoader_h
#define __CalibrationLoader_h

#include "MPXAlgo/MediPixAlgo.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"
#include "MAFTools/TimePixCalibrationHandler.h"

#define __calib_a_loaded   		0x01
#define __calib_b_loaded   		0x02
#define __calib_c_loaded   		0x04
#define __calib_t_loaded   		0x08
#define __calib_clock_loaded   	0x10

#define __calib_all_loaded 		0x1f

class CalibrationLoader {

public:

	CalibrationLoader(MediPixAlgo *);
	virtual ~CalibrationLoader() { };

	void SetCalibrationConfigFile_a(const char * f);
	void SetCalibrationConfigFile_b(const char * f);
	void SetCalibrationConfigFile_c(const char * f);
	void SetCalibrationConfigFile_t(const char * f);
	void SetCalibClk(double);

	void SetMaxFrameWidth(int mw) { m_maxWidth = mw; };
	void SetMaxFrameHeight(int mh) { m_maxHeight = mh; };

	void CalibrationInit();

	// Get the calibrated energy
	double CalculateAndGetCalibEnergy(int x, int y, int tot);
	double CalculateAndGetCalibEnergy(pair<int, int> pix, int tot);

	// Get the calibration handler
	MAFTools::TimePixCalibrationHandler * GetCalibrationHandler(){return m_calib;};
	double GetCalibClk(){return m_calibClock;};
	int GetMaxHeightCalib(){return m_maxHeight;};
	int GetMaxWidthCalib(){return m_maxWidth;};
	bool CalibIsOK(){return m_calibOK;};


private:

	MediPixAlgo * m_callerAlgo;

	string m_a_configFile;
	string m_b_configFile;
	string m_c_configFile;
	string m_t_configFile;
	int m_maxHeight;
	int m_maxWidth;
	double m_calibClock;

	// TOT calibration
	MAFTools::TimePixCalibrationHandler * m_calib;
	int  m_calibWord;
	bool m_calibOK;

};

#endif
